function MEEOsc = meanMEE2OscMEE(meeMean, param)
    
% convert to kepler 
kepmean = CoordConv.mee2coe(meeMean); 
      
% kepler to hill
hillmean = kep2hill(kepmean, param.mu);

% mean hill to osc hill
hillOsc = mean2osculatingHill(hillmean, param);

% hill to kepler
coeOsc = hill2kep(hillOsc, param.mu);

% kepler to MEE 
MEEOsc = CoordConv.kepler2MEOE([coeOsc(1), coeOsc(2), coeOsc(5),...
    coeOsc(4),coeOsc(3),coeOsc(6)]);


end

