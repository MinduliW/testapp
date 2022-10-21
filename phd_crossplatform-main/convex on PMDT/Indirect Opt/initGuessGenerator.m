function [L0exp,deltat,lamba_Eopt_linear] = initGuessGenerator(param)
%[x0, xf,xguess,controllimit] = initGuessGenerator(N, param)

[fsolveoptions, param.odeoptions] = fsolOpt(1e-10);
N = 1;
L0 = param.x0(6)+0; %linspace(0,-2*pi,16);

Lf = [0];

param.eclipse = false;

while ~(max(Lf)>param.xf(6)&& min(Lf)<param.xf(6))
    
    % set L0
    param.x0(6) = L0(N);
    
    % run l free energy optimal
    param.freeL = true;
    
    %zeros(1,6)'; %
    if N == 1
    lamba_x_t0 =zeros(1,6)';
    else
        lamba_x_t0= lamba_Eopt_linear';
    end
    
    lamba_Eopt_linear  =  fsolve(@(lambda0)propagateState(lambda0,@solve_x_linear, param),lamba_x_t0',fsolveoptions);
  
    [~,states] = ode45(@(t,x) solve_x_linear(t, x, param),[param.tvec],...
        [param.x0,param.m0,lamba_Eopt_linear],param.odeoptions);
    
    Lf(N) = states(end,6);
    
    if Lf(N)>param.xf(6)
        
        % go backwards in L0
        L0(N+1) =   L0(N) - 2*pi/30; 
    else
         % go forwards in L0
        L0(N+1) =   L0(N) + 2*pi/30; 
    end
    

        
        % If xf is not within the sequence, increase the sequence. 
        N = N+1;

     
 
end

% Interpolate through Lf to find apt L0.
L0exp = interp1((Lf), L0(1:end-1), (param.xf(6)));

% what is the change in time corresponding to this change in L0 ? 

eta_0 = sqrt(param.mu/param.x0(1)^3);


deltat = (L0exp - param.x0(6))/eta_0;

% testing
param.x0(6) = L0exp;



% run l free energy optimal
%param.freeL = true;


lamba_x_t0= lamba_Eopt_linear';


lamba_Eopt_linear  =  fsolve(@(lambda0)propagateState(lambda0,@solve_x_linear, param),lamba_x_t0',fsolveoptions);

[~,states] = ode45(@(t,x) solve_x_linear(t, x, param),[param.tvec],...
    [param.x0,param.m0,lamba_Eopt_linear],param.odeoptions);

if rem(param.xf(6)-states(end,6),pi)<1e-3
     param.xf(6)
    states(end,6)
    disp("Success : L of the debris reached from free L run");
    
end


