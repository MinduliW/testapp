function [dx,Ft,efficiency,L] = propagate2(t, x,cf, param)



le1 = cf(1);le2 = cf(2) ;lainc = cf(3);
leinc = cf(4);laRAAN = cf(5);lAOP = cf(6);


% Get A, B, C ,D and dB/dx for x.
[A,B,~, ~, ~] = getAandBMatrices(x,param);

coe = CoordConv.mee2coe(x);

%ArgofLat = wrapTo2Pi(coe(5) + coe(6));
    
[rr, vv] = CoordConv.po2pv(coe, param.mu);
   
hill = kep2hill(coe, param.mu);

% calculate mean hill
hillMean = osculating2meanHill(hill, param.mu, param.J2, param.Re);

coemean = hill2kep(hillMean, param.mu);



%% target tracking

[~, ~, xfKepOsc] = getPosition(param.t0+t*param.TU/86400, ...
       'H2AF15',3.986005e14);

xfKepOsc(1)= xfKepOsc(1)/param.LU;
hill = kep2hill(xfKepOsc, param.mu);
hillMean = osculating2meanHill(hill, param.mu, param.J2, param.Re);
xfKepOsc = hill2kep(hillMean, param.mu);
param.smaexp  = xfKepOsc(1);
param.RAANexp = xfKepOsc(4);
param.incexp = xfKepOsc(3);
param.eexp = xfKepOsc(2);

dL_doe = dVLaw(coemean, le1, le2, lainc, leinc, laRAAN, lAOP, param);

G  = gauss(coemean, param);
U = G*dL_doe;

Ft = -param.T/x(7)*U/norm(U);


L = lyapunov(coemean,cf,param);

dv= sqrt(L);

efficiency = -norm(Ft)*norm(U)/(2*dv);

% if efficiency >param.lcutoff
%     Ft = [0;0;0]; 
% end

if isnan(Ft)
    'here'
end


%% J2
F_J2 = J2acceleration(x, param)';


%% Drag 

%density = Density_nrlmsise00(coe(1)-param.Re, t + param.t0, param);

rho = param.dragfactor*Density_HP(coe(1) - param.Re);

v = sqrt(param.mu/coe(1));
a_d = 0.5*rho*param.Cd*param.Area*v^2/x(7);

% Drag acceleration is in the opposite direction to the velocity.
drag_accel_ECI = -a_d*vv/norm(vv);

% obtain the ECI to RTN matrix.
hhat = cross(rr, vv)/norm(cross(rr, vv));
rhat = rr/norm(rr);

% Now convert drag_accel_ECI to drag_accel_RTN
drag_accel_RTN = [rhat'; cross(hhat, rhat)' ; hhat']*drag_accel_ECI;

%drag_term = B*drag_accel_RTN;



if param.eclipses == true
     period =  2*pi*sqrt(coe(1)^3/param.mu);
    
    TA = wrapTo2Pi(coe(end));
    
    % on a grid of theta, determine the start and end (hence centre) of the
    % eclipse for the orbit.
    theta = linspace(TA,TA+2*pi,25);
    time = linspace(t,t+period,25);
    
    % for each theta, get the eclipse.
    for i = 1: length(theta)
        ecl(i) = isEclipse(time(i), coe(1),coe(2), coe(3), coe(4), coe(5), theta(i), param);
    end
    
    % get eclipse thetas
    eclipseindx = find(ecl ==0);
    
     if isempty(eclipseindx)
        theta_centre = 0;
    else
        
        theta_centre = (theta(eclipseindx(1))+theta(eclipseindx(end)))/2;
    end
    
    theta_centre_opp = wrapTo2Pi(theta_centre+pi);
    
     ArgofLat = wrapTo2Pi(coe(5) + coe(6));
    
    if ( ArgofLat>= theta_centre-pi/4) && (ArgofLat <= theta_centre+pi/4)
        eta = 0;
    elseif ( ArgofLat>= theta_centre_opp-pi/4) && (ArgofLat <= theta_centre_opp+pi/4)
        eta = 0;
    else
        eta = 1;
    end
else
    eta = 1;
end

% if param.dutyRatio ~= 1
% if (ArgofLat>pi/2 && ArgofLat<=3*pi/2)
%     eta = 0;
%     
% else
%     eta = 1;
% end
% end

if param.Thrust == false
    eta = 0;
    
end

F = eta*Ft + F_J2;


dx(1:6,1) = A+B*F;

if sum(Ft) ~= 0
    dx(7,1) = -eta*param.T/param.Isp/param.g0;
else
    dx(7,1) = 0;
end
