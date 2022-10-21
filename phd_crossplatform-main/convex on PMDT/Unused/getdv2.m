function [dvCart,accelCart] = getdv2(states,timevec, param)
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
    
 
    % Get A, B, C ,D and dB/dx for x.
    [~,B,~, ~, ~] = getAandBMatrices(states(i, 1:6)',param);
    
    
    % Define alpha_1
    alpha1 = -B'*states(i,8:13)'/norm(B'*states(i,8:13)');
    
    gamma1 = (norm(B'*states(i,8:13)'));
     
    rho = param.TR - gamma1;
    
    if rho > 0
        gamma1 = 0;
    else
        gamma1 = param.m0/states(i,7);
    end
    
    uu = alpha1*gamma1;

    posvel = CoordConv.ep2pv(states(i,1:6), param.mu);
    
    
    % convert to mean
  
    
    
    rr = posvel(1:3);
    vv = posvel(4:6);
    
    % gotta go from RTN to ECI for dvs.
    hhat = cross(rr, vv)/norm(cross(rr, vv));
    rhat = rr/norm(rr);
    
    RTN = [rhat'; cross(hhat, rhat)' ; hhat'];
    
    uuCart(:,i) = RTN\uu;
    accelCart(:,i) = uuCart(:,i)*param.Tmax/states(i,7);
    
    
    
     if ~isreal(uuCart)
        'here'
    end
    
end

for i = 1: length(timevec)-1
    
 ax = param.Tmax*uuCart(1,i:i+1)./states(i:i+1,7)';
 ay = param.Tmax*uuCart(2,i:i+1)./states(i:i+1,7)';
 az = param.Tmax*uuCart(3,i:i+1)./states(i:i+1,7)';
 dvCart(1,i) = trapz([timevec(i), timevec(i+1)], ax);
 dvCart(2,i) = trapz([timevec(i), timevec(i+1)], ay);
 dvCart(3,i) = trapz([timevec(i), timevec(i+1)], az);
 
 
%  dvCart(1,i) = param.Tmax*uuCart(1,i)*param.dt;
%  dvCart(2,i) = param.Tmax*uuCart(2,i)*param.dt;
%  dvCart(3,i) = param.Tmax*uuCart(3,i)*param.dt;
%  
%  

end
dvCart(:,i+1) = [0;0;0];
