function [x] = CoptSolve2(As,Bs,constants,xg,N,nlx, param)

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

[r, res] = mosekopt('symbcon');

prob.c   = [zeros(1,6*(N+1)+3*N), ones(1, N)];


% Dynamics and BCs 
Aeq3 = [];Beq3 = [];
for iters = 1: N
    
    ubd = 6*iters;
    lbd = ubd - 5;
    
    A = As(lbd:ubd,:);
    B = Bs(lbd:ubd, :);
    
    
    dynamicsArray = [zeros(6, 6*(iters-1)), A, -eye(6), zeros(6,6*(N-iters)),...
        zeros(6, 3*(iters-1)), B, zeros(6, 3*(N-iters)), zeros(6,N)];
    
    Aeq3 = [Aeq3; dynamicsArray];
    
    %zeros(1,6)'; %
    constantsarray =  -constants(lbd:ubd) + A*xg(iters,1:6)' + B*xg(iters,7:9)';
    
    Beq3 = [Beq3; constantsarray];
end

% initial position.
Aeq1 = ([eye(6), zeros(6,10*N)]);Beq1 = param.x0';

Aeq2 = ([zeros(6,6*N), eye(6), zeros(6,4*N)]); % ; % 
Beq2 = param.xf'; 
 
Aeq = [Aeq1; Aeq2; Aeq3];Beq = [Beq1; Beq2; Beq3];

prob.a = sparse(Aeq); % store matrix efficiently as a sparse matrix

% Lower and Upper Equality Constraints i.e. Ax=b
prob.blc = Beq;
prob.buc = Beq;


%compute spacecraft mass

if param.method == 1
    
    [~,mass] = ode45(@(t,m) massder(t,m,xg(:,7:9), param),[param.tvec],...
        [param.m0],odeset('RelTol', 1e-13,'AbsTol',1e-13));
    
    %         mdot = mass(i-1)*norm(xg(i-1,7:9)/param.Isp/param.g0);
    %         mass(i) = mass(i-1)-mdot*(param.dt);
    
    %figure; plot(mass);
    
else
    mass(1) = param.m0;
    for i = 2:N+1
        
        dv = norm(xg(i-1,7:9));
        mass(i) = mass(i-1)/exp(dv/param.Isp/param.g0);
    end
end



%% define a sphere around a given node in xg.
% eta = 0.01; 
% posradius = eta; 
% velradius = eta; 
% 
% xupper = [];
% xlower =[]; 
% for i = 1: N+1
%     xupper = [xupper, xg(i,1:3)+posradius,xg(i,4:6)+velradius];
%     xlower = [xlower, xg(i,1:3)-posradius,xg(i,4:6)-velradius];
% end
% 
% 
% if param.method == 1
%     uu= [];
%       uradius = param.Tmax./mass(1:end-1)'; %param.clim'; %
% 
%     for i = 1: N
%       
%       uu = [uu,uradius(i),uradius(i),uradius(i)];
%    
%     end
%     
%     prob.blx = [ xlower, -uu,  zeros(1, N)];
%     prob.bux = [xupper,  uu, uradius];
% 
% else
%     if param.eclipse == true
%         
%         for i = 1:N
%               v(i) = getEclipseCart(param.tvec(i), xg(i,1:6)', param);
%               
%               if v(i)>0
%                   v(i) = 1;
%               end
%               
%         end
%   
%         
%         dvmaglim = (1-v)*param.Tmax./mass(1:end-1)*param.dt;
%     
%     else
%          dvmaglim = param.Tmax./mass(1:end-1)*param.dt;
%     
%     end
%     
%    
%     dvmagarray =[];
%     for i = 1: N
%         dvmagarray = [dvmagarray, dvmaglim(i),dvmaglim(i),dvmaglim(i)];
%     end
%     
%     prob.blx = [ xlower, -dvmagarray,  zeros(1, N), zeros(1, 6*(N+2))];
%     prob.bux = [xupper,  dvmagarray,  dvmaglim, ones(1, 6*(N+2))];
%     
%     
% end


%% Setting blc and bux by the nonlinear index.
prob.blx  =[];
prob.bux =[];

[limits] = minvar(0.01)';

for i = 1: N
    prob.bux  = [prob.bux , xg(i,1:3)+limits(i,1:3),xg(i,4:6)+limits(i,4:6)];
    prob.blx = [prob.blx, xg(i,1:3)-limits(i,1:3),xg(i,4:6)-limits(i,4:6)];
end
prob.bux  = [prob.bux , xg(end,1:3),xg(end,4:6)];
prob.blx = [prob.blx, xg(end,1:3),xg(end,4:6)];


%% or to go with the dvlimit
%dvmaglim = param.Tmax./mass(1:end-1)*param.dt;

for i = 1:N
    v(i) = getEclipseCart(param.tvec(i), xg(i,1:6)', param);
    
    if v(i)>0
        v(i) = 1;
    end
    
end

dvmaglim = (1-v)*param.Tmax./mass(1:end-1)*param.dt;


dvmagarray =[];
for i = 1: N
    dvmagarray = [dvmagarray, dvmaglim(i),dvmaglim(i),dvmaglim(i)];
end

prob.blx = [ prob.blx, -dvmagarray,  zeros(1, N)];
prob.bux = [prob.bux ,  dvmagarray,  dvmaglim];



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
