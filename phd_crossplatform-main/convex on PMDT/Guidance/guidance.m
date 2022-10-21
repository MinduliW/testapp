function [ LTFinalPos,TOF,startmass] = guidance(a0, RAAN0_t0, inc0, ...
    incf,  omega0, e0, driftorbit,leg1,leg2, waitTimevec,param)
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
[times_guide, MEEguide,dv,fuelburnt, TOF] = propTrajwithWaitingLyp(a0, RAAN0_t0, inc0, ...
    incf,  omega0, e0, leg1, leg2,waitTimevec,driftorbit, param);

coemean_guide =[]; coe_guide = []; 

for j = 1: 1: length(times_guide)
    

    coec = CoordConv.mee2coe(MEEguide(j,1:6));
    
    coe_guide = [coe_guide; coec];
    
    % convert to hill
    hill = kep2hill(coec, param.mu);
    
    % calculate mean hill
    hillMean = osculating2meanHill(hill, param.mu, param.J2, param.Re);
    
    cmean = hill2kep(hillMean, param.mu);
    
    coemean_guide = [coemean_guide; cmean];
end


LTFinalPos = coe_guide(end,:);
startmass = MEEguide(end,7);

fprintf('Guidance \n')
fprintf('Fuel burnt (kg) : %f  \n', fuelburnt*param.MU);
fprintf('dv: %f m/s \n', dv*param.LU/param.TU);
fprintf('TOF (d) : %f  \n', TOF*param.TU/86400);%89.733129  
fprintf('--------------------------------- \n')