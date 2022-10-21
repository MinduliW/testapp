function [times, states,dv,fuelburnt,TOF] = propTrajwithWaitingLyp(a0, RAAN0_t0, inc0, ...
    incf,  omega0, e0, leg1, leg2,waitTimevec,x, param)
addpath('/Users/minduli/Astroscale_ADR/Main/Applications/Guidance/SlimLocoche');

%% Leg 1

mintheta = 0;
diffold = 100;
param.tracktarget =false;
param.followEdelbaum = true;
param.ef = 1e-4;

for t = 0:0.0001:2*pi
    diff = getx0(t, a0, RAAN0_t0, inc0,param);
    
    if abs(diff) < abs(diffold)
        mintheta = t;
    end
    
    diffold = diff;
end

x0 = CoordConv.kepler2MEOE([a0, e0,omega0,RAAN0_t0,inc0,mintheta]);


% calculate eclipse centres 
for i = 1: length(leg1.a)

 n = sqrt(param.mu/leg1.a(i)^3);
 TA = mintheta + n*leg1.t(i);    
coe = [leg1.a(i), e0, leg1.inc(i), omega0, leg1.RAAN, TA]; 
centre(i) = eclipse_centre(leg1.t(i),coe,param);

end

figure; plot(centre);
param.Thrust = true;
param.Edelsma = leg1.a;
param.Edelinc = leg1.inc;
param.EdelRAAN = leg1.RAAN;

if length(leg1.t) >1
    
    cf(1)= 0.1;%e1
    cf(3) = 1;%a,i
    cf(5) = 0.45;%a,RAAN
    
    
    cf(2) = 0.5;
    cf(4) = 1-cf(3);
    cf(6) = 1;
    [tl1,sl1] = ode45(@(t,x) propagate1(t,x,cf,leg1.t,centre,param),...
        [0, leg1.t(end)],[x0';param.m0]);
else
    tl1 = leg1.t(1);
    sl1 = [x0';param.m0]';
end


fprintf('Errors after Leg 1 \n')
errorcalc(tl1, sl1, param.Edelsma(end),param.Edelinc(end),...
    param.EdelRAAN(end),param);
fprintf('___________________________\n')

dv1 = param.Isp*param.g0*log(param.m0/sl1(end,end));


%% Waiting Orbit

param.followEdelbaum = false;

param.Edelsma = leg1.a(end)*ones(size(waitTimevec));
param.Edelinc = leg1.inc(end)*ones(size(waitTimevec));
param.EdelRAAN =param.EdelRAAN(end) -param.k*x(1)^7*cos(x(2))*(waitTimevec-leg1.t(end)); %rad
%param.ef = 1e-4;

if length(waitTimevec) >1
    
    param.Thrust = false;
    
    [tw,sw] = ode45(@(t,x) propagatew(t,x,param),[waitTimevec(1),...
        waitTimevec(end)],sl1(end,:));
    
else
    tw = waitTimevec(1);
    sw = sl1(end,:);
end


fprintf('Errors after waiting \n')

errorcalc(tw, sw, param.Edelsma(end),param.Edelinc(end),...
    param.EdelRAAN(end),param);
fprintf('___________________________\n')


%% Leg 2


% cf(1)= 0.1;%e1
% cf(3) = 1;%a,i
% cf(5) = 0.85;%a,RAAN
% 
% 
% cf(2) = 0.1;
% cf(4) = 1-cf(3);
% cf(6) = 1;

%coestartl2 = CoordConv.mee2coe(sw(end,:));

% calculate eclipse centres 
% for i = 1: length(leg2.a)
%     
%     n = sqrt(param.mu/leg2.a(i)^3);
%     TA = mintheta + n*leg2.t(i);
%     coe = [leg2.a(i), coestartl2(2), leg2.inc(i), coestartl2(5), leg2.RAAN, TA];
%     centre2(i) = eclipse_centre(leg2.t(i),coe,param);
%     
% end

param.Thrust = true;
% to do: make stop based on orbital plane achievement.

opt   = odeset('Events', @(t,y) stopCriterionTargetTracking(t,y,cf, param));
[tl2,sl2] = ode45(@(t,x) propagate2(t,x,cf,param),...
    [leg2.t(1), leg2.t(end)+15*24*60*60/param.TU],sw(end,:),opt);


for t = 1: length(tl2)
    [~,~,~, L(t)] = stopCriterionTargetTracking(tl2(t),...
        sl2(t,:),cf, param);
end

figure;
plot(sqrt(L));

%+10*24*60*60
dv2 = param.Isp*param.g0*log(sw(end,end)/sl2(end,end));

fprintf('Errors after Leg2 \n');
[~, ~, xfKepOsc] = getPosition(param.t0+tl2(end)*param.TU/86400, ...
       'H2AF15',3.986005e14);
xfKepOsc(1)= xfKepOsc(1)/param.LU;
hill = kep2hill(xfKepOsc, param.mu);
hillMean = osculating2meanHill(hill, param.mu, param.J2, param.Re);
xfKepOsc = hill2kep(hillMean, param.mu);
param.smaexp  = xfKepOsc(1);
param.RAANexp = xfKepOsc(4);
param.incexp = xfKepOsc(3);
param.eexp = xfKepOsc(2);

errorcalc(tl2,sl2,param.smaexp,param.incexp,...
   param.RAANexp ,param);
% 
% errorcalc(tl2,sl2, param.Edelsma(end),param.Edelinc(end),...
%     param.EdelRAAN(end),param);
fprintf('___________________________\n')


%% Totals

times = [tl1;tw; tl2];
states = [sl1;sw; sl2];

dv = dv1 + dv2;

fuelburnt = param.m0 - sl2(end,end);
TOF = times(end);

rmpath('/Users/minduli/Astroscale_ADR/Main/Applications/Guidance/SlimLocoche');

end
