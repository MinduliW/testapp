function [x,fval] = CoptSolveConeProg(As,Bs,constants,xguess,N,param)


% FUNCTION NAME:
%   CoptSolveConeProg
%
% DESCRIPTION:
%   Function runs coneprog on the given problem
%
% INPUT:
%   As - (double []) STM
%   Bs - (double []) Gamma
%   constants - (double []) Constant part of the state
%   xguess - (double []) guess for the trajectory and control
%   N - (double) number of nodes 
%   param - (struct) problem parameters
%
% OUTPUT:
%   x - (double []) Solution (trajectory+control)
%   fval - (double []) Objective function value

 %% Objective function
    
    f = [zeros(1,6*(N+1)+3*N), ones(1, N)]'; %*1/param.Omega^2;
    
    %% Thrust constraints
    
    b = zeros(3,1);
    gamma =0;
    for iters = 1: N
        
        A = [zeros(3, 6*(N+1)), zeros(3, 3*(iters-1)), eye(3),...
            zeros(3, 3*(N-iters)),zeros(3,N)];
        d = [zeros(6*(N+1)+3*N, 1); zeros(iters-1,1); ones(1, 1); zeros(N-iters,1)];

        socConstraints(iters) = secondordercone(A,b,d,gamma);
        
    end
    
    %% boundary constraints
    
    %compute spacecraft mass
mass(1) = param.m0;

for i = 2:N+1
    
    dv = norm(xguess(i-1,7:9));
    mass(i) = mass(i-1)/exp(dv/param.Isp/param.g0);
end


    dvmaglim = param.Tmax./mass(1:end-1)*param.dt;

dvmagarray =[];
for i = 1: N
    dvmagarray = [dvmagarray, dvmaglim(i),dvmaglim(i),dvmaglim(i)];
end


param.xmin = param.xlower; %[-xdist*ones(1,3), -xvel*ones(1,3)]; % %
param.xmax = param.xupper; % [xdist*ones(1,3), xvel*ones(1,3)]; %inf*ones(1,6); %;

param.xmin =[-1.0514   -0.5505   -1.1532   -0.8514   -0.4433   -0.9270];
param.xmax = [ 1.0478    0.5498    1.1627    0.8430    0.4434    0.9344];
lb = [ repmat(param.xmin, [1,N+1]), -dvmagarray,  zeros(1, N)];
ub = [repmat(param.xmax, [1,N+1]),  dvmagarray,  dvmaglim];


%     lb = [-inf*ones(1,(N+1)*6+3*N),  zeros(1, N)];
%     ub = [inf*ones(1,(N+1)*6+3*N),  inf*ones(1, N)];
%     
   
    %% Dynamics
    Aeq3 = [];
    Beq3 = [];
    for iters = 1: N
        
        ubd = 6*iters;
        lbd = ubd - 5;
        
        A = As(lbd:ubd,:);
        B = Bs(lbd:ubd, :);
        
        dynamicsArray = [zeros(6, 6*(iters-1)), A, -eye(6), zeros(6,6*(N-iters)),...
            zeros(6, 3*(iters-1)), B, zeros(6, 3*(N-iters)), zeros(6,N)];
        
        Aeq3 = [Aeq3; dynamicsArray];
        
        constantsarray = -constants(lbd:ubd) + A*xguess(iters,1:6)' + B*xguess(iters,7:9)';
        
        Beq3 = [Beq3; constantsarray];
    end
    
  
    
     %% Boundary condtions.
      % initial position.
    Aeq1 = ([eye(6), zeros(6,10*N)]);
    Beq1 = param.x0'; 
    
    Aeq2 = ([zeros(6,6*N), eye(6), zeros(6,4*N)]);
    Beq2 = param.xf';
    
      
    Aeq = [Aeq1; Aeq2; Aeq3];
    Beq = [Beq1; Beq2; Beq3];
    
    
    %% Solving
    
    options = optimoptions('coneprog','Display','iter');
    [x,fval] = coneprog(f,socConstraints,[],[],Aeq,Beq,lb,ub,options);
    
end
