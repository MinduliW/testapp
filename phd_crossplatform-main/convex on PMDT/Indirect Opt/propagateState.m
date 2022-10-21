function F = propagateState(lambda0,propFunction, param)

% FUNCTION NAME:
%   propagateState
%
% DESCRIPTION:
%   Propagates dynamics and calculates the residual of the final state.
%
% INPUT:
%   lambda0 - (double []) initial costates 
%   propFunction - (function []) dynamics function
%   param - (struct) Problem parameters 
%
% OUTPUT:
%   F - (double []) Residual
%

[~,states] = ode45(@(t,x) propFunction(t, x, param),...
    [param.t0, param.tf],[param.x0,param.m0, lambda0],param.odeoptions);

rev = floor(states(end,6)/(2*pi));

% convert states to mean.states(end,1:6); %
MEEMean= oscMEE2meanMEE(states(end,1:6), param);

MEEMean(6) = MEEMean(6)+ rev*2*pi;
% 
F = MEEMean-param.xf(1:6);

if param.freeL == true
    F = [MEEMean(1:5)-param.xf(1:5),states(end,end)];
end

end