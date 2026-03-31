function out = Dynamic_system_2D

out{1} = @init;
out{2} = @fuN_e_valval;
% out{3} = @jacobian;
% out{4} = @jacobianp;
% out{5} = @hessians;
% out{6} = @hessiansp;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function [tspan, X0, options] = init
handles = feval(Bending_system_2D);
X0=[0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 100];

% --------------------------------------------------------------------------

% Dynamical functions
function dXdt = fuN_e_valval(t, kmrgd, par_C0, par_zeta, par_kappa, par_tau_a)

    a_M = kmrgd(1);
    a_c = kmrgd(2);
    
    % Building blocks
    C = @(A) sqrt((1 - ((2.* a_c)./A - 1).^2))./(2*sqrt(a_c));
    C_dA = @(A) (2 .* a_c-A) ./ ...
        (2 * A.^2 .* sqrt(A - a_c));
    C_da_c = @(A) - 1 ./ (2 * A .* sqrt(A - a_c));
    
    da_Mdt = - 2 .* par_kappa .* ((C(a_M) - par_C0).^2 - (C((1 - a_M)) - par_C0).^2 + ...
        2 .* a_M .* (C(a_M) - par_C0) .* C_dA(a_M) - 2 .* (1 - a_M) .* (C((1 - a_M)) - par_C0) .* C_dA((1 - a_M)));
    
    da_cdt = - (4 * par_kappa .* (a_M .* (C(a_M) - par_C0) .* C_da_c(a_M) + ...
        (1 - a_M) .* (C((1 - a_M)) - par_C0) .* C_da_c((1 - a_M))) +  par_zeta*sqrt(pi)./sqrt(a_c));
    
    dXdt = [1/par_tau_a * da_Mdt; da_cdt];


% --------------------------------------------------------------------------
function jac = jacobian(t, kmrgd, par_C0, par_zeta, par_kappa, par_tau_a)
% --------------------------------------------------------------------------
function jacp = jacobianp(t, kmrgd, par_C0, par_zeta, par_kappa, par_tau_a)
% --------------------------------------------------------------------------
function hess = hessians(t, kmrgd, par_C0, par_zeta, par_kappa, par_tau_a)
% --------------------------------------------------------------------------
function hessp = hessiansp(t, kmrgd, par_C0, par_zeta, par_kappa, par_tau_a)
%---------------------------------------------------------------------------
function tens3  = der3(t, kmrgd, par_C0, par_zeta, par_kappa, par_tau_a)
%---------------------------------------------------------------------------
function tens4  = der4(t, kmrgd, par_C0, par_zeta, par_kappa, par_tau_a)
%---------------------------------------------------------------------------
function tens5  = der5(t, kmrgd, par_C0, par_zeta, par_kappa, par_tau_a)