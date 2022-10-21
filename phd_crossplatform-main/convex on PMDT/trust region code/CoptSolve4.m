function [x] = CoptSolve4(As,Bs,constants,xg,N,TR, param)

% FUNCTION NAME:
%   CoptSolve
%
% DESCRIPTION:
%   Function runs MOSEK on the given problem
%
% INPUT:
%   As - (double []) STM
%   Bs - (double []) Gamma
%   constants - (double []) Constant part of the state
%   xg - (double []) guess for the trajectory and control
%   N - (double) number of nodes 
%   param - (struct) problem parameters
%
% OUTPUT:
%   x - (double []) Solution (trajectory+control)


clear prob;

paramMosek=[];
paramMosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-14;
paramMosek.MSK_DPAR_INTPNT_TOL_PFEAS = 1e-14;
paramMosek.MSK_IPAR_INFEAS_REPORT_AUTO = 1;
paramMosek.MSK_IPAR_INFEAS_REPORT_LEVEL = 3;

xcoef = zeros(1,6*(N+1));
dvcoef = zeros(1,3*N);
dvmagcoef = ones(1, N);
scoef =  zeros(1,6*(N+1));
TRcoef = zeros(1,6*(N+1));
% s = x- xbar btw and TRcoef is the trust region variable.

prob.c   = [xcoef, dvcoef, dvmagcoef, scoef,TRcoef];

% Dynamics and BCs 
 
Aeq3 = [];Beq3 = [];
s_part = zeros(6, 6*(N+1)); % for Aeq1,Aeq2, Aeq3 .
TR_part = zeros(6, 6*(N+1));

for iters = 1: N
    
    ubd = 6*iters;
    lbd = ubd - 5;
    
    A = As(lbd:ubd,:);
    B = Bs(lbd:ubd, :);
    
   
    dynamicsArray = [zeros(6, 6*(iters-1)), A, -eye(6), zeros(6,6*(N-iters)),...
        zeros(6, 3*(iters-1)), B, zeros(6, 3*(N-iters)), zeros(6,N),s_part,TR_part];
    
    constantsarray =  -constants(lbd:ubd) + A*xg(iters,1:6)' + B*xg(iters,7:9)';
    
    Aeq3 = [Aeq3; dynamicsArray];
    Beq3 = [Beq3; constantsarray];
end

% initial position.
Aeq1 = ([eye(6), zeros(6,10*N),s_part,TR_part]);
Beq1 = param.x0';

Aeq2 = ([zeros(6,6*N), eye(6), zeros(6,4*N),s_part,TR_part]); % ; % 
Beq2 = param.xf'; 
 
% introduce s = x - xbar.
% in Ax = b form its xbar = x-s.
Aeq4 = [];Beq4 = [];
Aeq5 = []; Beq5 = [];
for iters = 1: N+1
    
    A = [zeros(6,6*(iters-1)), eye(6), zeros(6,6*(N+1-iters)), ...
        zeros(6,4*N),...
        zeros(6, 6*(iters-1)) , -eye(6), zeros(6,6*(N+1-iters)),TR_part];
    B = xg(iters,1:6)';
    
    Aeq4 = [Aeq4; A];
    Beq4 = [Beq4; B];
    
    A2 = [zeros(6, 10*N+6), zeros(6, 6*(N+1)),...
        zeros(6, 6*(iters-1)), eye(6), zeros(6, 6*(N+1-iters))];
    
    B2 = TR(1:6);
    
    Aeq5 = [Aeq5; A2];
    Beq5 = [Beq5; B2];
    
end


% introduce t = eta
Aeq = [Aeq1; Aeq2; Aeq3; Aeq4;Aeq5];
Beq = [Beq1; Beq2; Beq3; Beq4;Beq5];

prob.a = sparse(Aeq); % store matrix efficiently as a sparse matrix

% Lower and Upper Equality Constraints i.e. Ax=b
prob.blc = Beq;
prob.buc = Beq;


%compute spacecraft mass

%% eclipsing s
for i = 1:N
    if param.eclipse == true
          mee = CoordConv.vec2mee(xg(i,1:3), xg(i,4:6),param.mu);   
          v(i) = 1- getEclipse(param.tvec(i),mee, param);

       % v(i) =  1- getEclipseCart(param.tvec(i), xg(i,1:6)', param);
    else
        v(i)=1;
    end
    
    
end

%%
if param.method == 1
    
    [~,mass] = ode45(@(t,m) massder(t,m,xg(:,7:9), param),[param.tvec],...
        [param.m0],odeset('RelTol', 1e-13,'AbsTol',1e-13));
    

else
    mass(1) = param.m0;
    for i = 2:N+1
        
        dv = norm(xg(i-1,7:9));  
        
        mass(i) = mass(i-1)/exp(dv/param.Isp/param.g0);


    end
end

%% Setting blc and bux by the nonlinear index.
prob.blx  =[];
prob.bux =[];

TR = TR';

for i = 1: N+1
    prob.bux  = [prob.bux , xg(i,1:3)+TR(1:3),xg(i,4:6)+TR(4:6)];
    prob.blx = [prob.blx, xg(i,1:3)-TR(1:3),xg(i,4:6)-TR(4:6)];
end

%dvmaglim = v*param.Tmax./mass(1:end-1)*param.dt;

% find dv limit by integrating Tmax/M using trapezium rule.
dvmaglim =[];
for i = 1:N
    a = param.Tmax/mass(i);
    b = param.Tmax/mass(i+1);
    h = param.dt; 
    area = (a+b)/2*h;
    dvmaglim(i) = v(i)*area;

end


dvmagarray =[];
for i = 1: N
    dvmagarray = [dvmagarray, dvmaglim(i),dvmaglim(i),dvmaglim(i)];
end

prob.blx = [ prob.blx, -dvmagarray,  zeros(1, N)];
prob.bux = [prob.bux ,  dvmagarray,  dvmaglim];

% add the limits for s and t too.
prob.blx  = [prob.blx , -Beq5', Beq5'];
prob.bux  = [prob.bux ,  Beq5', Beq5'];



%% 

prob.cones.type = zeros(N +6*(N+1),1);

l = 1:3:3*N;
m = 2:3:3*N;
n = 3:3:3*N;

prob.cones.sub = [];

prob.cones.subptr =[];
index = 6*N+6;
for i = 1:N % For each node
  
    prob.cones.sub  = [prob.cones.sub 9*N+6+i index+l(i) index+m(i) index+n(i)]; % Add the cone constraint
end
prob.cones.subptr = 1:4:4*(N); % Specify beginning index of each cone


% introduce cone contraints on s. stating |s| < t
for i = 1: 6*(N+1)
    
    tindex = 6*N+6 + 4*N+6*N+6 +i;
    sindex = 6*N+6 + 4*N+i;
    
     prob.cones.sub  = [prob.cones.sub tindex sindex]; % Add the cone constraint

end

startindx = 4*(N)+1;
endindx = startindx+2*6*(N+1)-2;

prob.cones.subptr = [prob.cones.subptr, startindx:2:endindx];





[r,res]=mosekopt('minimize',prob,paramMosek);

% Display the primal solution.

x= res.sol.itr.xx';


end
