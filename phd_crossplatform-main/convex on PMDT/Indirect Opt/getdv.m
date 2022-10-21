function [dvCart,accelCart] = getdv(states,timevec, param)
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
    
    uu= tau*alpha; 
    
 
    posvel = CoordConv.ep2pv(states(i,1:6), param.mu);
      
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
  
 dvCart(1,i) = param.Tmax*uuCart(1,i)*param.dt;
 dvCart(2,i) = param.Tmax*uuCart(2,i)*param.dt;
 dvCart(3,i) = param.Tmax*uuCart(3,i)*param.dt;
 

end



dvCart(:,i+1) = [0;0;0];
