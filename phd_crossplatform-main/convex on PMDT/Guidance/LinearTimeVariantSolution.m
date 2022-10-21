function lamba_x_ts =LinearTimeVariantSolution(param)

% FUNCTION NAME:
%   LinearTimeVariantSolution
%
% DESCRIPTION:
%   Calculates linear analytical energy optimal costates
%
% INPUT:
%   param - (struct) Problem parameters 
%
% OUTPUT:
%   lamba_x_ts - (double []) Costates
%

param.ts = 0;
options =param.odeoptions;
% Solve equation 22 to get phi.
[~,phis] = ode45(@(t,phi) solve_phi(t,phi, param),...
    [param.ts param.tf],eye(12),options);

phi_tf_ts = reshape(phis(end, :), [12,12]);


phi_11 = phi_tf_ts(1:6, 1:6); 
phi_12 = phi_tf_ts(1:6,7:12);

% Now calculate z. 
[~,z] = ode45(@(t,z) solve_z(t,z, param),...
    [param.ts param.tf],zeros(1,6),options);


z_tf_ts = z(end, :)';


% WrapT02pi on L
value = (phi_11*param.x0' + z_tf_ts);
%value(6) = wrapTo2Pi(value(6));

% calculate lambda_ts 
lamba_x_ts = (phi_12)\(param.xf' - value);


function dphi = solve_phi(t, phi, param)

% The shape of phi and dphi is an issue here. 
% reshape phi
phi = reshape(phi, [12,12]);

% get x from t. 
x = param.x0 + (param.xf - param.x0)/(param.tf - param.ts)*(t - param.ts);

% Generate new A B C D matricies.
[~,B,C, ~] = getAandBMatrices(x',param); 

% Construct the F matrix at tf. 
F_t = [C -param.Tmax/param.m0*(B*B'); zeros(6,6) -C'];


% Calculate dphi.
dphi = F_t*phi;

% Reshape dphi.
dphi = reshape(dphi, [144,1]);



function dz = solve_z(t,z, param)

x = param.x0 + (param.xf - param.x0)/(param.tf - param.ts)*(t - param.ts);


% Generate new A B C D matricies.
[~,~,C, D] = getAandBMatrices(x',param);


dz = C*z + D;

