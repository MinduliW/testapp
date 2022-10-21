clear,clc, close all;
addpath('Edelbaum');
addpath('Guidance');
addpath('Indirect Opt');
addpath('Dynamics');
addpath('Utility');

% check the platform cause my code is now cross platform/annoying 

if ismac
    param.mac = true;
elseif ispc
    param.mac = false;
end

if param.mac == true
    addpath('/Users/minduli/mosek/9.3/toolbox/r2015a');
    addpath('/Users/minduli/libraries_mice/mice/src/mice/');
    addpath('/Users/minduli/libraries_mice/mice/lib');
    addpath('SGP4routines_NAIF')
    cspice_furnsh('SGP4routines_NAIF/kernelmac.txt')
else

    addpath('H:\Libraries\Mosek\10.0\toolbox\r2017a');
    addpath('H:/Libraries/libraries_mice/mice/src/mice/');
    addpath('H:/Libraries/libraries_mice/mice/lib')
    cspice_furnsh('SGP4routines_NAIF/kernelwin.txt')
end



problemParam; 
param.Topt = true;param.addDrag = false;

%% RAAN MATCHING METHOD 

[driftorbit, Total_dv1,TOF,Omegadiff ,leg1,leg2, waitTimevec, targetRAAN, wait,~]...
    =RAANmatchingmethod(v0,vf, inc0, incf, param);

% set x0 and xf. 



%% 

% number of segments of convex 
convSeg = 500; 

startinx = 1; 
for i = 1: convSeg

    % divide the time into segments 
    endindx = round(i*length(leg1.t)/convSeg); 

    x0cartMEAN = [leg1.a(startinx), 0, leg1.inc(startinx), leg1.RAAN(startinx), 0,0 ]; 
    xfcartMEAN = [leg1.a(endindx), 0, leg1.inc(endindx), leg1.RAAN(endindx), 0,sqrt(param.mu/leg1.a(endindx)^3)*leg1.t(endindx) ]; 

    [N, param] =SCvxParam(x0cartMEAN,xfcartMEAN, leg1.t(startinx:endindx), param);

    xguess =[]; 

    N = 0; 
    for j = startinx: endindx
        N = N+1;
        theta(j) = sqrt(param.mu/leg1.a(j)^3)*leg1.t(j);
        coord =  [leg1.a(j), 0, leg1.inc(j), leg1.RAAN(j), 0, theta(j)];

        arglat = wrapTo2Pi(theta(j));


        % calculate components of dv
        beta = leg1.beta(j);

        if (arglat >= pi/2 && arglat < 3*pi/2)
            beta = abs(beta);
        else
            beta = -abs(beta);
        end
        
        dv = leg1.dV(j)*[0, cos(beta), sin(beta)];
        xguess(j,:) = [CoordConv.po2pv(coord,param.mu)',dv];
        

    end

    fileID = fopen('iter.txt','w');
    fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n',xguess');
    fclose(fileID);


    % ah gotta get dat dv components
    param.method = 2; 
    param.eclipse = false; 
    param.nu = 0.01;

    % write problem parameters to give to cpp.
    writeparams(param,N); 


    [dynamics,control,eta,fval] = runSCvx(xguess, N, param);
    plottrustregion(param)


end



figure;
plot3(rs(:,1), rs(:,2), rs(:,3), 'MarkerSize',1);
%% SCVX

