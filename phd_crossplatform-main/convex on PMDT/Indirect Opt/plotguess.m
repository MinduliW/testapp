function plotguess(lambda0, propfunction, accelfunction,tvec, param)

% FUNCTION NAME:
%   plotguess
%
% DESCRIPTION:
%   Plots the initial guess 
%
% INPUT:
%   lambda0 - (double []) initial costates 
%   propfunction - (function []) dynamics function
%   accelfunction - (function []) acceleration function
%   tvec - (double []) time vector
%   param - (struct) Problem parameters 
%

% Use ODE45 to propagate the states in time.
solve_xfunction = @(t,x) propfunction(t, x, param);

options = param.odeoptions;
[timevec,states] = ode45(solve_xfunction,...
    [tvec],[param.x0,param.m0,lambda0],options);

% find the size of states
[~, cols] = size(states); 

% Plot initial point 
rs = CoordConv.ep2pv(param.x0(1:6),param.mu);
rs = rs(1:3);

% Plot destination point.
rf =CoordConv.ep2pv(param.xf(1:6),param.mu);
rf = rf(1:3);

for i = 1: length(timevec)
    r = CoordConv.ep2pv(states(i,1:6), param.mu);
    positions(i,:) = r(1:3)';
    velocities(i,:)= r(4:6)';
    
end


% Position vs Time plot
figure;  hold on;
plot3(rs(1), rs(2), rs(3), '*');
plot3(rf(1), rf(2), rf(3), '*');
p = plot3(positions(:,1),positions(:,2),positions(:,3), '-');

%axis equal;
plot_latex(p, 'x (LU)', 'y (LU)','z (LU)', 'Initial Guess' ...
    ,{'$x_s$', '$x_f$',  'Initial Guess Trajectory' });

% 
% figure; hold on;
% p = plot3(velocities(:,1),velocities(:,2),velocities(:,3), '-');
% 
% plot_latex(p, 'vx (LU/TU)', 'vy(LU/TU)','vz(LU/TU)', 'Initial Guess' ...
%     ,{ });



%% Acceleration and Thrust

% Get the accelerations and plot them.
[uu,taus] = accelfunction(states, timevec, param);


figure;
subplot(3,1,1); hold on;
p = plot(timevec, taus);
plot_latex(p, 'time (TU)', '','', 'Thrust magnitude vs Time' ...
    ,{ '$\tau$'});

subplot(3,1,2);
p = plot(timevec, uu);
plot_latex(p, 'time (days)', 'Acceleration','', 'Acceleration vs Time' ...
    ,{'$a_r$', '$a_t$', '$a_n$'});

subplot(3,1,3);
p = plot(timevec*param.TU, states(:,7));

plot_latex(p, 'time (s)', 'Mass (m0)','', 'Mass vs Time' ...
    ,{});

