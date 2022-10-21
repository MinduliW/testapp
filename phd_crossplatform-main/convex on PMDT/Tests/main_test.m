clear,clc, close all;
addpath('/Users/minduli/Astroscale_ADR/Main/Core');
addpath('/Users/minduli/mosek/9.3/toolbox/r2015a');

constants; 
problem;

%% RAAN MATCHING METHOD 

% Initial guesses
param.Topt = true;param.addDrag = true;
x0 = [v0-0.1, inc0-0.001,0.01];

lb = [min([v0,vf])-1000, min([inc0,incf])-20*pi/180, 0 ];
ub = [max([v0,vf]), max([inc0,incf])+20*pi/180, inf];

if param.Topt == false
    ub(end) = 50*param.TU; % kinda need to limit this to be reasonable othewise the guidance will take eons to run
end

[x,~] = fmincon( @(x) fitnessfunction(x, param),x0,[],[],[],[],lb,ub,...
    @(x) condition2(x,param) , optimoptions('fmincon','Display','none'));

[Total_dv1,TOF,Omegadiff ,leg1,leg2, waitTimevec, targetRAAN, wait,fuelburnt]   = ...
    outcomes(x,0 , param);

fprintf('Edelbaum \n')
fprintf('Fuel burnt: %f kg \n', fuelburnt);
fprintf('dv: %f m/s \n', Total_dv1);
fprintf('TOF (d) : %f  \n', TOF/param.TU);
fprintf('--------------------------------- \n')

%% GUIDANCE - Ruggeiro 
% param.eps = 1;
% param.etae = 0; % param.dutyRatio;
% param.etaa = 0; %param.dutyRatio;
% param.etai = 0; %param.dutyRatio;
% param.etaRAAN =0; % param.dutyRatio;
% param.guidance = true;
% 
% % [times_guide, MEEguide,dv,fuelburnt,TOF] = propTrajwithWaiting(a0, RAAN0_t0,...
% %     inc0, omega0, e0, incf, leg1, leg2,waitTimevec,x, param);
% 
% 

%% Guidance - Locoche
[times_guide, MEEguide,dv,fuelburnt,~, TOF] = propTrajwithWaitingLyp(a0, RAAN0_t0, inc0, ...
    incf,  omega0, e0, leg1, leg2,waitTimevec,x, param);

fprintf('Errors with Guidance \n')
[coemean_guide,coe_guide]  = errorcalc(times_guide, MEEguide, af,incf,leg2.RAAN(end),param);

coemean1 = coemean_guide;
times1 = times_guide;

%plottrajwaiting(a0, af, inc0, x, TOF ,leg1,leg2,coemean_guide, times_guide,wait,param)

fprintf('Guidance \n')
fprintf('Fuel burnt (kg) : %f  \n', fuelburnt);
fprintf('delta v : %f  \n', dv);
fprintf('dv: %f m/s \n', dv);
fprintf('TOF (d) : %f  \n', TOF);%89.733129  
fprintf('--------------------------------- \n')

% need to reduce this error as much as possible.