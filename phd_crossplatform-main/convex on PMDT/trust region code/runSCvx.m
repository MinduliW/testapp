function [dynamics,control,eta,fval] = runSCvx(xguess, N, param)

% FUNCTION NAME:
%   runSCvx
%
% DESCRIPTION:
%   Runs the SCVx algorithm to obtain successive convexification outcomes
%
% INPUT:
%   N - (double) number of nodes 
%   xguess - (double []) initial guess of the trajectory
%   param - (struct) Problem parameters 
%


% Plot initial point
rs = CoordConv.ep2pv(param.x0(1:6),param.mu);
rs = rs(1:3);

% Plot destination point.
rf =CoordConv.ep2pv(param.xf(1:6),param.mu);
rf = rf(1:3);

figure;hold on;
plot3(rs(1), rs(2), rs(3), '*');
plot3(rf(1), rf(2), rf(3), '*');

param.x0(1:6) = CoordConv.ep2pv(param.x0(1:6),param.mu);
param.xf(1:6) =CoordConv.ep2pv(param.xf(1:6),param.mu);

p = plot3(xguess(:,1), xguess(:,2),xguess(:,3));


legnd{1} = 'Departure';
legnd{2} = 'Target';
legnd{3} = 'Initial guess';
distance = 100;
diter = 0;
while diter <20 && distance >0.05
    
    diter = diter +1;
    
    legnd{diter+3} = strcat('Iteration:', num2str(diter));
    !./run
    
    %%
    
    [As, Bs,constants,TR] =extractfromtext('outputcpp.txt','trustRegion.txt');
    
    [x] = CoptSolve4(As,Bs,constants,xguess,N,TR,param);
    
    %[x,fval] = CoptSolveConeProg(As,Bs,constants,xguess,N,param);
    
    dynamics = reshape(x(1:6*(N+1)), [6, N+1]);
    control = reshape( x(6*(N+1)+1: 6*(N+1)+3*N), [3, N]);
    control = [control'; zeros(3,1)'];
    eta = x(6*(N+1)+3*N+1:6*(N+1)+3*N+N);
    
    fval(diter) = sum(eta);
    
    xnew = [dynamics',control];
  
        
    p = plot3(dynamics(1,:), dynamics(2,:),dynamics(3,:), '+-');
    
    
    % write xguess into a text file.
    fileID = fopen('iter.txt','w');
    fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n',xnew');
    fclose(fileID);
    
    %% Find difference between current and guess
    
    xguesspath = xguess(:,1:3);
   
    xnewpath = xnew(:,1:3);
    
    etaguess = 0;
    for val = 1: length(xguesspath)
         etaguess =  etaguess+ norm(xguess(val,7:9));
        distance(val) = norm(xguesspath(val,:) - xnewpath(val,:));
    end
    
    
    distance = sum(distance);
    


    distance

    %distance = norm(xguesspath(end,1:3)- param.xf')
    xguess = xnew;
    
    if param.method == 1
        fprintf("total u = %f \n", sum(eta)*param.LU/param.TU/param.TU) ;

    else
        
        fprintf("total dv = %f \n", sum(eta)*param.LU/param.TU) ;
        
    end
    
    
end

plot_latex(p, 'x(m)', 'y(m)','z (m)', '' ,legnd)


figure;hold on;
p = plot3(dynamics(1,:), dynamics(2,:),dynamics(3,:),'+-');
plot3(rs(1), rs(2), rs(3), '*');
plot3(rf(1), rf(2), rf(3), '*');

plot_latex(p, 'x(m)', 'y(m)','z (m)', 'Final solution',{})

if param.method == 2
    quiver3(dynamics(1,:),dynamics(2,:),dynamics(3,:),control(:,1)', control(:,2)', control(:,3)',1, 'LineWidth', 2); % plot the delta-Vs
end



end
