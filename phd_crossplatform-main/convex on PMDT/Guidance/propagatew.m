function [dx] = propagatew(t, x, param)

% Get A, B, C ,D and dB/dx for x.
[A,B,~, ~, ~] = getAandBMatrices(x,param);


%% J2
F_J2 = J2acceleration(x, param)';


F = F_J2;


dx(1:6,1) = A+B*F;
dx(7,1) = 0;
  

