clear,clc, close all;
addpath('Edelbaum');
addpath('Guidance');
addpath('Indirect Opt');
addpath('Dynamics');
addpath('Utility')
addpath('/Users/minduli/mosek/9.3/toolbox/r2015a');
addpath('/Users/minduli/libraries_mice/mice/src/mice/')
addpath('/Users/minduli/libraries_mice/mice/lib')
cspice_furnsh('/Users/minduli/Astroscale_ADR/Main/SGP4routines_NAIF/kernel.txt')

problemParam; 
param.Topt = true;param.addDrag = false;

%% RAAN MATCHING METHOD 

[driftorbit, Total_dv1,~,Omegadiff ,leg1,leg2, waitTimevec, targetRAAN, wait,~]...
    =RAANmatchingmethod(v0,vf, inc0, incf, param);

%% GUIDANCE - Ruggeiro 

[LTFinalPos,TOF,mass] = guidance(a0, RAAN0_t0, inc0, ...
    incf,  omega0, e0, driftorbit,leg1,leg2, waitTimevec,param);

%% Convex Optimization 

%% SCVX problem setup 

proxHours = 3;


[N,param] = SCvxParam(LTFinalPos,param.t0+TOF*param.TU/86400,proxHours, param);
% 
%% Scaling for COPT

param.addJ2 = true;
% acceleration or dv as control?
param.method = 2; %acceleration = 1, %dv = 2; 
param.eclipse = true;
param.cs =1000;
param.ct =1;
param.nu = 0.5;

% write problem parameters to give to cpp.
writeparams(param,N); 

% Initial guess 
 [L0exp,param.deltat0,lamba_Eopt_linear] = initGuessGenerator(param);
 
param.x0(6) = L0exp;
param.freeL = false;
[xguess,param.clim] = initGuess(N,lamba_Eopt_linear, param);

%%

[dynamics,control,eta,fval] = runSCvx(xguess, N, param);
plottrustregion(param)

 %% forward propagate the obtained solution

rs = CoordConv.ep2pv(param.x0(1:6),param.mu)';
rf = CoordConv.ep2pv(param.xf(1:6),param.mu)';

figure;  hold on;
plot3(rs(1), rs(2), rs(3), '*');
plot3(rf(1), rf(2), rf(3), '*');

tvec = linspace(param.t0, param.tf,N+1);

if param.method ==2
    for i = 1:N
        
        % add dv to the current state
        rs(4:6) = rs(4:6)+control(i,:);
        
        % store the state with added dv
        statesAfterdv(i,:) = rs;
        
        % propagate from current node to the next.
        [t,propstate] = ode45(@(t,x) propagateCart(t, x, param),...
            [tvec(i),tvec(i+1)],[rs], odeset('RelTol', 1e-13,'AbsTol',1e-13));
        
        
        % make rs the propagated state.
        rs =propstate(end,:);
        
    end
    
    statesAfterdv(i+1,:) = rs;
else
    param.control = control;
    param.time = tvec;

    [t,statesAfterdv] = ode45(@(t,x) propagateCart(t, x, param),...
        [param.tvec],[rs],odeset('RelTol', 1e-13,'AbsTol',1e-13));
    
    
end

p = plot3(statesAfterdv(:,1),statesAfterdv(:,2),statesAfterdv(:,3), 'x-');

if param.method == 1
    plot_latex(p, 'x (LU)', 'y (LU)','z (LU)', 'Final solution forward validation' ...
        ,{'$x_s$', '$x_f$',  'trajectory'});
    
    figure;
    p = plot(param.tvec*param.TU/(60*60), [eta,0]*param.LU/param.TU/param.TU);
    plot_latex(p, 'time (h)', '$u (m/s^2)$ ','', 'Final solution' ...
        ,{});
    
    figure; hold on;
    p = plot(param.tvec*param.TU/(60*60), [control(:,1)']*param.LU/param.TU/param.TU);
    p = plot(param.tvec*param.TU/(60*60), [control(:,2)']*param.LU/param.TU/param.TU);
    p = plot(param.tvec*param.TU/(60*60), [control(:,3)']*param.LU/param.TU/param.TU);
    plot_latex(p, 'time (h)', '$u (m/s^2)$','', 'Final solution' ...
        ,{'$u_x$', '$u_y$', '$u_z$'});
    
    
    
    [~,mass] = ode45(@(t,m) massder(t,m,control, param),[param.tvec],...
        [param.m0],odeset('RelTol', 1e-13,'AbsTol',1e-13));
    deltav = param.Isp*param.g0*log(param.m0/mass(end))*param.LU/param.TU;
    fprintf("total dv = %f \n", deltav) ;
    
    
    
else

    quiver3(statesAfterdv(:,1),statesAfterdv(:,2),statesAfterdv(:,3),control(:,1), control(:,2), control(:,3), 0.5, 'LineWidth', 2); % plot the delta-Vs
    plot_latex(p, 'x (LU)', 'y (LU)','z (LU)', 'Final solution forward validation' ...
        ,{'$x_s$', '$x_f$',  'trajectory' ,'$\Delta v$'});
    
    % calculate total dv.
    deltav =  sum(eta)*param.LU/param.TU;
    fprintf("total dv = %f \n",deltav) ;
    
    for i = 1:N
        
        % convert to MEOE
        mee = CoordConv.vec2mee( dynamics(1:3,i), dynamics(4:6,i),param.mu);   
        v(i) = getEclipse(param.tvec(i),mee, param);

    end
    
        
    figure; subplot(2,1,1); hold on;
      p = plot(param.tvec*param.TU/(60*60), [(1-v),0]);
    p = plot(param.tvec*param.TU/(60*60), [eta,0]*param.LU/param.TU);
    plot_latex(p, 'time (h)', '$\Delta v$ (m/s)','', 'Final solution' ...
        ,{'$\nu$ (Eclipse)', '$\eta$ (Control Mag.)'});
    
    subplot(2,1,2); hold on;
    p = plot(param.tvec*param.TU/(60*60), [control(:,1)']*param.LU/param.TU);
    p = plot(param.tvec*param.TU/(60*60), [control(:,2)']*param.LU/param.TU);
    p = plot(param.tvec*param.TU/(60*60), [control(:,3)']*param.LU/param.TU);
    plot_latex(p, 'time (h)', '$\Delta v$ (m/s)','', 'Final solution' ...
        ,{'$\Delta v_x$', '$\Delta v_y$', '$\Delta v_z$'});
end

% calculate fuel expenditure
mf = param.m0/exp( sum(eta)/param.Isp/param.g0);
fuelmass = (param.m0-mf)*param.MU;

% convert dv to RTN
for i = 1: N
    
    rr = dynamics(1:3,i);
    vv = dynamics(4:6,i);
    
    hhat = cross(rr, vv)/norm(cross(rr, vv));
    rhat = rr/norm(rr);
    
    RTN = [rhat'; cross(hhat, rhat)' ; hhat'];
    
    dvRTN(i,:) = RTN*control(i,:)';
    
    
end

figure; hold on;
    p = plot(param.tvec(1:end-1)*param.TU/(60*60), [dvRTN(:,1)']*param.LU/param.TU);
    p = plot(param.tvec(1:end-1)*param.TU/(60*60), [dvRTN(:,2)']*param.LU/param.TU);
    p = plot(param.tvec(1:end-1)*param.TU/(60*60), [dvRTN(:,3)']*param.LU/param.TU);
    plot_latex(p, 'time (h)', '$\Delta v$ (m/s)','', '$\Delta v$ in RTN' ...
        ,{'$\Delta v_R$', '$\Delta v_T$', '$\Delta v_N$'});
 