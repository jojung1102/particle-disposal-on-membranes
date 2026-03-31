function out = Bending_system_2D

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
function dXdt = fuN_e_valval(t, kmrgd, par_C0, par_zeta, par_kappa)

    a_M = kmrgd(1);
    a_c = kmrgd(2);
    
    % Building blocks
    C = @(a) sqrt((1 - ((2.* a_c)./a - 1).^2))./(2*sqrt(a_c));
    C_da = @(a) (2 .* a_c-a) ./ ...
        (2 * a.^2 .* sqrt(a - a_c));
    C_da_c = @(a) - 1 ./ (2 * a .* sqrt(a - a_c));
    
    da_Mdt = - 2 .* par_kappa .* ((C(a_M) - par_C0).^2 - (C((1 - a_M)) - par_C0).^2 + ...
        2 .* a_M .* (C(a_M) - par_C0) .* C_da(a_M) - 2 .* (1 - a_M) .* (C((1 - a_M)) - par_C0) .* C_da((1 - a_M)));
    
    da_cdt = - (4 * par_kappa .* (a_M .* (C(a_M) - par_C0) .* C_da_c(a_M) + ...
        (1 - a_M) .* (C((1 - a_M)) - par_C0) .* C_da_c((1 - a_M))) +  par_zeta*sqrt(pi)./sqrt(a_c));
    
    dXdt = [da_Mdt; da_cdt];


% --------------------------------------------------------------------------
function jac = jacobian(t, kmrgd, par_C0, par_zeta, par_kappa)
% --------------------------------------------------------------------------
function jacp = jacobianp(t, kmrgd, par_C0, par_zeta, par_kappa)
% --------------------------------------------------------------------------
function hess = hessians(t, kmrgd, par_C0, par_zeta, par_kappa)
% --------------------------------------------------------------------------
function hessp = hessiansp(t, kmrgd, par_C0, par_zeta, par_kappa)
%---------------------------------------------------------------------------
function tens3  = der3(t, kmrgd, par_C0, par_zeta, par_kappa)
%---------------------------------------------------------------------------
function tens4  = der4(t, kmrgd, par_C0, par_zeta, par_kappa)
%---------------------------------------------------------------------------
function tens5  = der5(t, kmrgd, par_C0, par_zeta, par_kappa)