function centre = eclipse_centre(t,coe,param)
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
    centre = 0;
else
    
    centre = (theta(eclipseindx(1))+theta(eclipseindx(end)))/2;
end

end
