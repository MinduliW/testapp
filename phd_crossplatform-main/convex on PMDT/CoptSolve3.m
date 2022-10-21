function [x] = CoptSolve3(As,Bs,constants,xg,N,TR, param)

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

prob.c   = [zeros(1,6*(N+1)+3*N), ones(1, N),ones(1,6*N), ones(1,2*6)]; %+...

% Dynamics and BCs 
Aeq3 = [];Beq3 = [];
for iters = 1: N
    
    ubd = 6*iters;
    lbd = ubd - 5;
    
    A = As(lbd:ubd,:);
    B = Bs(lbd:ubd, :);
    
    dynamicsArray = [zeros(6, 6*(iters-1)), A, -eye(6), zeros(6,6*(N-iters)),...
        zeros(6, 3*(iters-1)), B, zeros(6, 3*(N-iters)), zeros(6,N), ...
        zeros(6, 6*(iters-1)), -eye(6), zeros(6, 6*(N-iters)) , zeros(6,6),zeros(6,6)];
    
    constantsarray =  -constants(lbd:ubd) + A*xg(iters,1:6)' + B*xg(iters,7:9)';
    
    Aeq3 = [Aeq3; dynamicsArray];
    Beq3 = [Beq3; constantsarray];
end

% initial position.
Aeq1 = ([eye(6), zeros(6,10*N), zeros(6, 6*N), eye(6), zeros(6,6)]);
Beq1 = param.x0';

Aeq2 = ([zeros(6,6*N), eye(6), zeros(6,4*N),zeros(6, 6*N), zeros(6,6), eye(6)]); % ; % 
Beq2 = param.xf'; 
 
Aeq = [Aeq1; Aeq2; Aeq3];Beq = [Beq1; Beq2; Beq3];

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

for i = 1: N
    prob.bux  = [prob.bux , xg(i,1:3)+TR(1:3,i)',xg(i,4:6)+TR(4:6,i)'];
    prob.blx = [prob.blx, xg(i,1:3)-TR(1:3,i)',xg(i,4:6)-TR(4:6,i)'];
end


prob.bux  = [prob.bux , xg(end,1:3)+TR(1:3,end)',xg(end,4:6)+TR(4:6,end)'];
prob.blx = [prob.blx, xg(end,1:3)-TR(1:3,end)',xg(end,4:6)-TR(4:6,end)'];


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


% figure;
% plot(dvmaglim)


m1 = 1;

if m1 == 1
dvmagarray =[];
for i = 1: N
    dvmagarray = [dvmagarray, dvmaglim(i),dvmaglim(i),dvmaglim(i)];
end

prob.blx = [ prob.blx, -dvmagarray,  zeros(1, N)];
prob.bux = [prob.bux ,  dvmagarray,  dvmaglim];


else
    
    for i = 1: N
        
        % upper limits.
        if norm(xg(i,7:9)+TR(7:9)) <  dvmaglim(i)
            prob.bux = [prob.bux , xg(i,7:9)+TR(7:9)];
            udv(i) = norm(xg(i,7:9)+TR(7:9));
            
        else
            prob.bux = [prob.bux ,  dvmaglim(i), dvmaglim(i), dvmaglim(i)];
            udv(i) = norm( dvmaglim(i));
        end
        
        if norm(xg(i,7:9)-TR(7:9)) <  dvmaglim(i)
            prob.blx = [prob.blx , xg(i,7:9)-TR(7:9)];
            
            if norm(xg(i,7:9)-TR(7:9))>0
                ldv(i) = norm(xg(i,7:9)-TR(7:9));
            else
                ldv(i) = 0;
            end
            
        else
            prob.blx = [prob.blx ,  -dvmaglim(i), -dvmaglim(i), -dvmaglim(i)];
            ldv(i) = 0;
        end
        
    end
    
    

    prob.bux  = [prob.bux ,udv];
    prob.blx = [prob.blx,ldv];
    
end





prob.bux  = [prob.bux ,ones(1, 6*(N+2))];
prob.blx = [prob.blx,zeros(1, 6*(N+2))];

%% 

prob.cones.type = zeros(N,1);

l = 1:3:3*N;
m = 2:3:3*N;
n = 3:3:3*N;

prob.cones.sub = [];
index = 6*N+6;
for i = 1:N % For each node
  
    prob.cones.sub  = [prob.cones.sub 9*N+6+i index+l(i) index+m(i) index+n(i)]; % Add the cone constraint
end
prob.cones.subptr = 1:4:4*(N); % Specify beginning index of each cone

[r,res]=mosekopt('minimize',prob,paramMosek);

% Display the primal solution.

x= res.sol.itr.xx';


end

