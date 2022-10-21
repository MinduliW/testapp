function [xguess,controllimit] = initGuess(N,lamba_Eopt_linear, param)


% FUNCTION NAME:
%   initGuess
%
% DESCRIPTION:
%   Calculates an initial guess for the RDV problem using the linear energy
%   optimal solution.
%
% INPUT:
%   N - (double) No of nodes
%
% OUTPUT:
%   xguess - (double []) Initial guess trajectory and control
%   controllimit - (double) Upper limit for the control magnitude
%


param.eclipse = false;

[fsolveoptions, param.odeoptions] = fsolOpt(1e-13);

%lamba_x_t0 =LinearTimeVariantSolution(param);

lamba_Eopt_linear  =  fsolve(@(lambda0)propagateState(lambda0,@solve_x_linear, param),lamba_Eopt_linear,fsolveoptions);

plotguess(lamba_Eopt_linear, @solve_x_linear, @accel_x_linear,param.tvec, param);

[timevec,states] = ode45(@(t,x) solve_x_linear(t, x, param),[param.tvec],...
[param.x0,param.m0,lamba_Eopt_linear],param.odeoptions);
states(end,6)
[dv,aaCart] = getdv(states,timevec, param);

for i = 1:length(timevec)
    
    if param.method == 1
        u(i,:)  = aaCart(:,i);
        normu(i) = norm(aaCart(:,i));
    else
        u(i,:)  = dv(:,i);
    end
    
    
    % convert states to cartesian
    statesC(i,1:6) = CoordConv.ep2pv(states(i,1:6), param.mu);
    
    %u(i,:) = [0,0,0];%REMOVE THIS IF CONTROL IN GUESS
    
    xguess(i,:) = [statesC(i,1:6), u(i,:)];
end

controllimit = param.Tmax./states(1:end-1,7);


fileID = fopen('iter.txt','w');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n',xguess');
fclose(fileID);


% Plot initial point
rs = CoordConv.ep2pv(param.x0(1:6),param.mu)';

rf =CoordConv.ep2pv(param.xf(1:6),param.mu)';

% Position vs Time plot
figure;  hold on;
plot3(rs(1), rs(2), rs(3), '*');
plot3(rf(1), rf(2), rf(3), '*');


if param.method == 1
    
    param.time = timevec;
    param.control = u;
    
    [t,propstate] = ode45(@(t,x) propagateCart(t, x, param),...
        [param.tvec],[rs],param.odeoptions);
    
    p = plot3(propstate(:,1),propstate(:,2),propstate(:,3), 'x-');
 
    plot_latex(p, 'x (LU)', 'y (LU)','z (LU)', 'Initial guess forward validation' ...
        ,{'$x_s$', '$x_f$',  'trajectory'});
    
    figure;
    p = plot(param.tvec*param.TU/(60*60), [normu]*param.LU/param.TU/param.TU);
    plot_latex(p, 'time (h)', '$u(m/s^2)$ ','', 'Initial Guess' ...
        ,{});
    
    
    
else
    
    for i = 1:N
        
        % add dv to the current state
        rs(4:6) = rs(4:6)+u(i,:);
        
        dvs(i) = norm(u(i,:));
        
        % store the state with added dv
        statesAfterdv(i,:) = rs;
        
        %update mass
        
        % propagate from current node to the next.
        [t,propstate] = ode45(@(t,x) propagateCart(t, x, param),...
            [param.tvec(i),param.tvec(i+1)],[rs],param.odeoptions);
        
        % make rs the propagated state.
        rs =propstate(end,:);
        
        
        statesAfterdv(i+1,:) = rs;
        
    end
    
    p = plot3(statesAfterdv(:,1),statesAfterdv(:,2),statesAfterdv(:,3), 'x-');
    quiver3(statesAfterdv(:,1),statesAfterdv(:,2),statesAfterdv(:,3),u(:,1), u(:,2), u(:,3), 0.6, 'LineWidth', 2); % plot the delta-Vs
    
    plot_latex(p, 'x (LU)', 'y (LU)','z (LU)', 'Initial guess forward validation' ...
        ,{'$x_s$', '$x_f$',  'trajectory' ,'$\Delta v$'});
    
    fprintf("total dv guess = %f \n", sum(dvs)*param.LU/param.TU) ;
    
    
    figure;
    p = plot(param.tvec*param.TU/(60*60), [dvs,0]*param.LU/param.TU);
    plot_latex(p, 'time (h)', '$\Delta v$ (m/s)','', 'Initial Guess' ...
        ,{});
    
end

end



% posradius = [ 0.198    0.106    0.2185]; %*(abs(sin(pi/(0.5*(N))*(NVec-4))));
% velradius = [ 0.16    0.083    0.173]; %*(abs(sin(pi/(0.5*(N))*(NVec-4))));
% 
% 
% xupper = [];
% xlower =[]; 
% for i = 1: N+1
%     xupper = [xupper, xguess(i,1:3)+posradius,xguess(i,4:6)+velradius];
%     xlower = [xlower, xguess(i,1:3)-posradius,xguess(i,4:6)-velradius];
% end
% 
