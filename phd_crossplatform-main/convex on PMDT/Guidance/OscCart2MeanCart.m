function posvelmean = OscCart2MeanCart(posvel, param)
    
% convert to kepler 
x = CoordConv.vec2orbElem(posvel(1:3),posvel(4:6), param.mu);
      
% kepler to hill
hill = kep2hill(x, param.mu);

% osc hill to mean hill
hillMean = osculating2meanHill(hill, param.mu, param.J2, param.Re);

% hill to kepler
coemean = hill2kep(hillMean, param.mu);

% kepler to cart 
[rrmean, vvmean] = CoordConv.po2pv(coemean, param.mu);

posvelmean = [rrmean; vvmean];

