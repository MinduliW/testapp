Tmax =1;

% Max Ts that work. 
Isp = 300;
g0 = 9.80665;
mu_earth = 3.986004418e14;
m0 = 100;

Rs = 695800e3;
Re = 6378.1363e3;


LU = Re; %1000e3; %
VU = sqrt(mu_earth/Re);
TU = LU/VU;
MU = m0;

param.LU = LU;
param.TU = TU;
param.MU = MU;


% Scale and store variables
param.Re = Re/LU;
param.RS = Rs/LU;

param.J2 =  1.08262668e-3; 
param.Isp = Isp/param.TU;
param.g0 = g0/(LU)*param.TU*param.TU;
param.Tmax = Tmax/MU/(LU)*TU^2;
param.m0 = m0/MU;
param.mu = mu_earth/LU^3*TU^2;
param.Rs = Rs/LU;

%% Positions and times

[xs, xf, TOF] = setCase('F');
param.ephemData = load('Debris.txt');
TOF = 1.0*24*60*60;
idd = 115;
ida = 82;
epochd = 23780;
RAN2 = []; RAAN1 = [];

start = 2.378651e4-TOF/2/(24*60*60);
endt = 2.378651e4+TOF/2/(24*60*60);
[xs(1:3), xs(4:6), EE1] = DebrisEphem(start, idd, param.ephemData,3.986004418e5, 1.08262668e-3, 6378.137);

[xf(1:3), xf(4:6), EE2] = DebrisEphem(endt, ida, param.ephemData, 3.986004418e5, 1.08262668e-3, 6378.137);

% convert from km to m.
xs = xs*1e3;
xf = xf*1e3;

time = linspace(start,endt,10);
for i = 1 : length(time)
    [~,~, EE1] = DebrisEphem(time(i), idd, param.ephemData,3.986004418e5, 1.08262668e-3, 6378.137);
    RAAN1(i) = EE1(4);
    [~,~, EE2] = DebrisEphem(time(i), ida, param.ephemData, 3.986004418e5, 1.08262668e-3, 6378.137);
    RAAN2(i) = EE2(4);
end

figure; hold on;
plot(time,RAAN1);
p = plot(time,RAAN2);
plot_latex(p, 'time', '$\Omega$','', '' ,{'Initial', 'Final'})

% % Scaling + CoordConv
param.x0 = CoordConv.vec2mee(xs(1:3)/LU,xs(4:6)/LU*TU,param.mu);
param.xf =  CoordConv.vec2mee(xf(1:3)/LU,xf(4:6)/LU*TU,param.mu);

% Time vector. Time is in TU
param.t0 =0; param.tf = param.t0 + TOF/TU;

period1 = 2*pi*sqrt(param.x0(1)^3/param.mu);
period2 = 2*pi*sqrt(param.xf(1)^3/param.mu);
period = 0.5*(period1+period2);
revolutions = round(param.tf/period);
param.xf(6) = param.xf(6) + (revolutions-1)*2*pi;


if revolutions > 1
    N = round(param.tf/period*25);
else
    N = 108;
end
