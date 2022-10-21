function [uu,taus,rho] = accel_x_linear(states,timevec, param)
% FUNCTION NAME:
%   accel_x_linear
%
% DESCRIPTION:
%   Function gives the acceleration and the dv of the trajectory
%
% INPUT:
%   states - (double []) x values 
%   timevec - (double []) time vector
%   param - (struct) problem parameters
%
% OUTPUT:
%   uu - (double []) acceleration vector
%   taus - (double []) thrust acceleration ratio
%   rho - (double []) switching function value 
%   dvCart - (double []) dv vector


for i = 1: length(timevec)
    
   % alpha_e = param.Tmax/states(i,7); 
   x = param.x0 + (param.xf - param.x0)/(param.tf - param.t0)*(timevec(i) - param.t0);

    
    % Get A, B, C ,D and dB/dx for x.
    [~,B,~, ~, ~] = getAandBMatrices(x',param);
    
    % Define alpha_1
    
    Bl = B'*states(i, 8:13)';
    alpha = -Bl/(norm(Bl));
    
    tau = norm(Bl);
    rho(i) = nan;
    
    taus(i) = tau;
    uu(:,i) = -Bl;
    
   
end




end
