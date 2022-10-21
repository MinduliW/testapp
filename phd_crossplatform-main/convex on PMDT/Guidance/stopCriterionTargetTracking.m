function [position,isterminal,direction,L] = stopCriterionTargetTracking(t,y, cf,param)

% convert final position to keplerian.
coe_end = CoordConv.mee2coe(y); 
hill = kep2hill(coe_end, param.mu);
hillMean = osculating2meanHill(hill, param.mu, param.J2, param.Re);
coemean_end = hill2kep(hillMean, param.mu);

[~, ~, xfKepOsc] = getPosition(param.t0+t*param.TU/86400,'H2AF15',3.986005e14);
xfKepOsc(1)= xfKepOsc(1)/param.LU;
hill = kep2hill(xfKepOsc, param.mu);
hillMean = osculating2meanHill(hill, param.mu, param.J2, param.Re);
xfKepOsc = hill2kep(hillMean, param.mu);
param.smaexp  = xfKepOsc(1);
param.RAANexp = xfKepOsc(4);
param.incexp = xfKepOsc(3);
param.eexp = xfKepOsc(2);

RAANgap =   abs(coemean_end(4)-param.RAANexp);
incgap = abs(coemean_end(3)-param.incexp);


plane_gap = (RAANgap+incgap)*180/pi;


L = lyapunov(coemean_end,cf,param);

dv = sqrt(L);
val = ~(dv<0.5);

% 
% 
% if L > 1e12
% %     'here'
% %     L
% end
% 
% dV = sqrt(L);



position = val; % The value that we want to be zero
isterminal = 1;  % Halt integration
direction = 0;   % The zero can be approached from either direction
end