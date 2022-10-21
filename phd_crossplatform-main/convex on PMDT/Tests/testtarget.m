clear, clc;

p2; 

period = 2*pi*sqrt(af^3/param.mu);

t = linspace(0,period,100);

rrs =[];
for i = 1: length(t)
[rr, vv, xfKepOsc] = getPosition(param.t0+t(i)*param.TU/86400, ...
    'H2AF15',mu, J2, Re);

rrs = [rrs,rr];
end

figure;
plot3(rrs(1,:),rrs(2,:),rrs(3,:))