function [mindist, maxdist] = trajlimits(param)

figure; hold on;

tol = 1e-13;
% propagate the initial orbit.
rs = CoordConv.ep2pv(param.x0(1:6),param.mu);

period = 2*pi*sqrt(param.mu/rs(1)^3);

[t,propstate1] = ode45(@(t,x) propagateCart(t, x, param),...
    [0, param.tf],[rs], odeset('RelTol', tol,'AbsTol',tol));

% find the minimum ECI coordinate

    
plot3(propstate1(:,1),propstate1(:,2), propstate1(:,3));


% Plot destination point.
rf =CoordConv.ep2pv(param.xf(1:6),param.mu);
period = 2*pi*sqrt(param.mu/rf(1)^3);

[t,propstate2] = ode45(@(t,x) propagateCart(t, x, param),...
    [0,  param.tf],[rf], odeset('RelTol', tol,'AbsTol',tol));

% find the max ECI coordinate
propstate = [propstate1;propstate2];


mindist= min(propstate(:,1:3));
maxdist = max(propstate(:,1:3));

plot3(propstate2(:,1),propstate2(:,2), propstate2(:,3));

end
