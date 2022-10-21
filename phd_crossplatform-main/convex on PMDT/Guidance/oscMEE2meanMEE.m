function MEEMean= oscMEE2meanMEE(MEEOsc, param)
    
% convert to kepler 
kepOsc = CoordConv.mee2coe(MEEOsc); 
      
% kepler to hill
 hillOsc = kep2hill(kepOsc, param.mu);

% mean hill to osc hill
hillmean = osculating2meanHill(hillOsc, param.mu, param.J2, param.Re);

% hill to kepler
coemean = hill2kep(hillmean, param.mu);

% kepler to MEE 
MEEMean = CoordConv.kepler2MEOE([coemean(1), coemean(2), coemean(5),...
    coemean(4),coemean(3),coemean(6)]);


end
