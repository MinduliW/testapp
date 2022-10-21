function initGuess2(N, param)

addpath('/Users/minduli/phd-thesis/Improving Indirect Methods/Final Matlab/fopt');


param.ts = param.t0;
param.xs = param.x0;

[fsolveoptions, param.odeoptions,param.bvpoptions] = fsolOpt(1e-13);


lamba_fsolve_linear =LinearTimeVariantSolution(param)';

if param.eclipse == true
    eclipseinit = true;
else
    eclipseinit = false;
end

param.eclipse = false; 

param.massvariation = false; 
lamba_fsolve_nonlinear  = fsolve(@(lambda0)propagateState(lambda0,@solve_x_eopt,...
    param),[lamba_fsolve_linear],fsolveoptions);

param.massvariation = true; 
lamba_fsolve_nonlinear  = fsolve(@(lambda0)propagateState(lambda0,@solve_x_eopt,...
    param),[lamba_fsolve_nonlinear],fsolveoptions);


%% Calculating Tau threshold. 

% Set k and taumax.
param.tau_max = 1; 

param.TR =get_Thrust_threshold(lamba_fsolve_nonlinear,@solve_x_eopt,@accel_eopt, param);


[time_energyopt_nonlin,states_energyopt_nl] = ode45(@(t,x) ...
    solve_x_eopt(t, x, param),[param.ts, param.tf],...
    [param.xs,param.m0, lamba_fsolve_nonlinear],param.odeoptions);

solinit.x = time_energyopt_nonlin';
solinit.y =  [states_energyopt_nl]';



%% Solving the fuel optimal problem.c
%k = 0.6666:0.1111:0.9999;
param.k =  0; % [0.333:0.333:0.999]
kmax = 0.99;

method =1; 

param.massvariation = true;
initiallamba = lamba_fsolve_nonlinear;
param.full = false; param.cm = false;
i = 1;
while  param.k  <= 1
    
   % param.k = k(i);
    disp(param.k);
    
    if i == 1 || method ~= 1 
        mainfunction = @(t,x) solveXFOpt(t, x, param);
        boundaryconditions = @(ya,yb) conditions(ya,yb,param);
        
        fuelOptimalSolution = bvp5c(mainfunction,boundaryconditions,solinit,param.bvpoptions);
        initiallamba = fuelOptimalSolution.y(8:end,1)';
        
        solinit.x = fuelOptimalSolution.x;
        solinit.y = fuelOptimalSolution.y;
        
    end
    
        
    if method == 1
    initiallamba  = fsolve(@(lambda0)propagateState(lambda0,@solveXFOpt,...
        param),initiallamba,fsolveoptions);
    end
    
    
    param.k = param.k +kmax/5;

    i = i+1;

end
toc

param.k = kmax;
param.full = true;
initiallamba  = fsolve(@(lambda0)propagateState(lambda0,@solveXFOpt,...
    param),initiallamba,fsolveoptions);



%% Introduce Eclipse.

initiallamba_e = initiallamba; 
if eclipseinit == true
    param.eclipse = true;
    
    epsilon = 0:0.1:1;
    
    for i = 1:length(epsilon)
        
        param.eps = epsilon(i);
        disp(epsilon(i));
        
        initiallamba_e = fsolve(@(lambda0)propagateState(lambda0,...
            @solveXFOpt,param),[initiallamba_e],fsolveoptions);
        
    end
end


%%
propForward(initiallamba_e',param)
 
[~,states] = ode45(@(t,x)solveXFOpt(t, x, param),[param.tvec],...
    [param.xs,param.m0, initiallamba_e],param.odeoptions);

endmass = states(end,7);
deltav = param.Isp*param.g0*log(param.m0/endmass)*param.LU/param.TU;

fprintf("total dv from indirect  = %f \n", deltav) ;
% 
% % find the difference between eopt trajectory and final trajectory.
% plotguess(initiallamba, @solveXFOpt, @accel_fopt,param.tvec, param);
% 
% [dv,aaCart] = getdv2(states,timevec, param);
% 
% 
% for i = 1:length(param.tvec)
%     
%    if param.method == 1
%         u(i,:)  = aaCart(:,i);
%         normu(i) = norm(aaCart(:,i));
%     else
%         u(i,:)  = dv(:,i);
%     end
%     
%     % convert states to cartesian
%     statesC(i,1:6) = CoordConv.ep2pv(states(i,1:6), param.mu);
%     
%     %u(i,:) = [0,0,0];
%     xguess(i,:) = [statesC(i,1:6), u(i,:)];
% end
% 
% 
% controllimit = param.Tmax./states(1:end-1,7);
% 
% figure;
% plot(controllimit)
% 
% 
% param.time = timevec;
% param.control = u;
% 
% fileID = fopen('iter.txt','w');
% fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n',xguess');
% fclose(fileID);
% 
% 
% % Plot initial point
% rs = CoordConv.ep2pv(param.x0(1:6),param.mu)';
% 
% rf =CoordConv.ep2pv(param.xf(1:6),param.mu)';
% 
% % Position vs Time plot
% figure;  hold on;
% plot3(rs(1), rs(2), rs(3), '*');
% plot3(rf(1), rf(2), rf(3), '*');
% 
% 
% if param.method == 1
%     
%     
%     [t,propstate] = ode45(@(t,x) propagateCart(t, x, param),...
%         [param.tvec],[rs],param.odeoptions);
%     
%     p = plot3(propstate(:,1),propstate(:,2),propstate(:,3), 'x-');
%  
%     plot_latex(p, 'x (LU)', 'y (LU)','z (LU)', 'Initial guess forward validation' ...
%         ,{'$x_s$', '$x_f$',  'trajectory'});
%     
%     figure;
%     p = plot(param.tvec*param.TU/(60*60), [normu]*param.LU/param.TU/param.TU);
%     plot_latex(p, 'time (h)', '$u(m/s^2)$ ','', 'Initial Guess' ...
%         ,{});
%     
%  
% else
%     
%     for i = 1:N
%         
%         % add dv to the current state
%         rs(4:6) = rs(4:6)+u(i,:);
%         
%         dvs(i) = norm(u(i,:));
%         
%         % store the state with added dv
%         statesAfterdv(i,:) = rs;
%         
%         %update mass
%         
%         % propagate from current node to the next.
%         [t,propstate] = ode45(@(t,x) propagateCart(t, x, param),...
%             [param.tvec(i),param.tvec(i+1)],[rs],param.odeoptions);
%         
%         % make rs the propagated state.
%         rs =propstate(end,:);
%         
%         
%         statesAfterdv(i+1,:) = rs;
%         
%     end
%     
%     p = plot3(statesAfterdv(:,1),statesAfterdv(:,2),statesAfterdv(:,3), 'x-');
%     quiver3(statesAfterdv(:,1),statesAfterdv(:,2),statesAfterdv(:,3),u(:,1), u(:,2), u(:,3), 0.6, 'LineWidth', 2); % plot the delta-Vs
%     
%     plot_latex(p, 'x (LU)', 'y (LU)','z (LU)', 'Initial guess forward validation' ...
%         ,{'$x_s$', '$x_f$',  'trajectory' ,'$\Delta v$'});
%     
%    %fprintf("total dv guess = %f \n", sum(dvs)*param.LU/param.TU) ;
%     
%     
%     figure;
%     p = plot(param.tvec*param.TU/(60*60), [dvs,0]*param.LU/param.TU);
%     plot_latex(p, 'time (h)', '$\Delta v$ (m/s)','', 'Initial Guess' ...
%         ,{});
%     
% end
end

