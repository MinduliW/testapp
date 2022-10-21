function xdot = propagateCart(t,x,param)
 
% FUNCTION NAME:
%   propagateCart
%
% DESCRIPTION:
%   Propagates dynamics, USed for testing only.
%
% INPUT:
%   t - (double) time
%   x - (double []) x in Cartesian coordinates 
%   param - (struct) Problem parameters 
%
% OUTPUT:
%   xdot - (double []) time derivative of state x. 
%


r = sqrt(x(1)^2+ x(2)^2 + x(3)^2);

xdot(1) = x(4); 
xdot(2) = x(5); 
xdot(3) = x(6); 
xdot(4) = -param.mu*x(1)/r^3; 
xdot(5) = -param.mu*x(2)/r^3; 
xdot(6) = -param.mu*x(3)/r^3; 

if param.addJ2 == true
    
    xdot(4) =xdot(4) -param.mu*x(1)/r^3*1.5*param.J2*(param.Re/r)^2*(1- 5*x(3)^2/r^2);
    xdot(5) =xdot(5)-param.mu*x(2)/r^3*1.5*param.J2*(param.Re/r)^2*(1- 5*x(3)^2/r^2);
    xdot(6) =xdot(6) -param.mu*x(3)/r^3*1.5*param.J2*(param.Re/r)^2*(3- 5*x(3)^2/r^2);
 


end

if param.method == 1
    
    ux = interp1(param.time, param.control(:,1), t);
    uy = interp1(param.time, param.control(:,2), t);
    uz = interp1(param.time, param.control(:,3), t);
    
    
    xdot(4) =xdot(4)+ux;
    xdot(5) =xdot(5)+uy;
    xdot(6) =xdot(6)+uz;
    
    
end


xdot = xdot';
end
