function [xguess,controllimit] = initGuess3(N, param)

addpath('/Users/minduli/phd-thesis/Low thrust Di Wu/Working version');
tol = 1e-13;
fsolveoptions = optimoptions('fsolve','Display','iter');
fsolveoptions.FunctionTolerance = tol; 
fsolveoptions.FunctionTolerance = tol;
fsolveoptions.MaxIterations = 5e3;

param.beta = param.Tmax/(param.Isp*param.g0);
param.odeoptions =  odeset('RelTol', tol,'AbsTol',tol);
param.bvpoptions = bvpset('Stats','on','AbsTol', tol, 'RelTol', tol);

%% Linearised outcome 

% In this section equations 22, 23 and 26 are used to obtain the linear
% result of the low thrust trajectory. 
lamba_0_linearised  =  [LinearTimeVariantSolution(param)',0];
fprintf('%s %.4f %.4f %.4f %.4f %.4f %.4f %.4f  \n' , 'lambda_s from analytical eq : ' , lamba_0_linearised)

%% Runing Algorithm 1 

% Define initial Epsilon
param.eps = 1;

% Set target epsilon
eps_tar = 0; 

%% Plotting the results.

% Define ts
param.ts = param.tf -param.eps*(param.tf - param.t0);



fprintf('%s %.4f \n', 'Epsilon = ', param.eps);
[lamba_fsolve,fval,exitflag,output]  = fsolve...
        (@(lambda0) propagateState1(lambda0, param),lamba_0_linearised,fsolveoptions);
   
disp(lamba_fsolve);

%plotresults(lamba_fsolve, param);

lambas= [lamba_fsolve]; 

tic
% Repeat till param.eps = eps_tar 
while param.eps ~= eps_tar 
    
    
     if exitflag <= 0 
        %d_epsilon = 0.5*d_epsilon;
        disp('Error: Initial solution did not work')
    else 
        d_epsilon = min([0.05, param.eps-eps_tar]); 
    end
   
    % calculate new epsilon
    param.eps = param.eps - d_epsilon;

    param.ts = param.tf -param.eps*(param.tf - param.t0);

    fprintf('%s %.4f \n', 'Epsilon = ', param.eps);
    % use fsolve to generate the new solution.
    [lamba_fsolve,fval,exitflag,output]  = fsolve...
        (@(lambda0) propagateState1(lambda0, param),lamba_fsolve,fsolveoptions);
    
    disp(lamba_fsolve);

    lambas = [lambas;lamba_fsolve];
    
    
end

toc 

plotresults(lamba_fsolve, param);
toc

xguess = 0;
controllimit = 0;
