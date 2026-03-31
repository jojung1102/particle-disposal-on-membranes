function out = Bending_system_4D

out{1} = @init;
out{2} = @fuNeval;
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
function [tspan, X0,options] = init
handles = feval(Bending_system_4D);
X0=[0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 100];

% --------------------------------------------------------------------------

% Dynamical functions
function dXdt = fuNeval(t, kmrgd, par_C0, par_c, par_zeta, par_kappa, par_n_f, par_n_d, par_epsilon, par_t_n, par_t_A)

    a_M = kmrgd(1);
    a_c = kmrgd(2);
    n_fM = kmrgd(3);
    n_dM = kmrgd(4);

    gamma = 1 ./ (1 + exp(-par_epsilon));
    
    C = @(A) sqrt((1 - ((2.* a_c)./A - 1).^2))./(2*sqrt(a_c));
    C_dA = @(A) - (A - 2 .* a_c) ./ ...
        (2 * A.^2 .* sqrt(A - a_c));
    C_da_c = @(A) - 1 ./ (2 * A .* sqrt(A - a_c));
    
    n_M = n_dM + n_fM;
    n_V = par_n_f + par_n_d - n_M;
    C_M = par_c * n_M./a_M;
    C_V = par_c * n_V./(1 - a_M);
    
    entropic_contribution = log((a_M .* ((1 - a_M) - n_V)) ./ ((a_M - n_M) .* (1 - a_M)));
    
    % Change in variables
    da_Mdt = entropic_contribution - 2 .* par_kappa .* ((C(a_M) - par_C0 - C_M).^2 - (C((1 - a_M)) - par_C0 - C_V).^2 + ...
        2 .* a_M .* (C(a_M) - par_C0 - C_M) .* (C_dA(a_M)+C_M/a_M) - 2 .* (1 - a_M) .* (C((1 - a_M)) - par_C0 - C_V) .* (C_dA(1 - a_M)+C_V/(1 - a_M)));
    
    da_cdt = - (4 * par_kappa .* (a_M .* (C(a_M) - par_C0 - C_M) .* C_da_c(a_M) + ...
        (1 - a_M) .* (C((1 - a_M)) - par_C0 - C_V) .* C_da_c((1 - a_M))) +  par_zeta*sqrt(pi)./sqrt(a_c));

    dn_fMdt = -log((gamma .* n_fM .* (1 - a_M-n_V)) ./ ((par_n_f - n_fM) .* (a_M-n_M))) + 4 * par_kappa * par_c * ((C(a_M) - par_C0 - C_M) - (C(1 - a_M) - par_C0 - C_V));
    dn_dMdt = -log((n_dM .* (1 - a_M-n_V)) ./ ((par_n_d - n_dM) .* (a_M-n_M))) + 4 * par_kappa * par_c * ((C(a_M) - par_C0 - C_M) - (C(1 - a_M) - par_C0 - C_V));

    dXdt = [da_Mdt/par_t_A; da_cdt; dn_fMdt/par_t_n; dn_dMdt/par_t_n]; 

% --------------------------------------------------------------------------
function jac = jacobian(t, kmrgd, par_C0, par_c, par_zeta, par_kappa, par_n_f, par_n_d, par_epsilon, par_t_n, par_t_A)
% --------------------------------------------------------------------------
function jacp = jacobianp(t, kmrgd, par_C0, par_c, par_zeta, par_kappa, par_n_f, par_n_d, par_epsilon, par_t_n, par_t_A)
% --------------------------------------------------------------------------
function hess = hessians(t, kmrgd, par_C0, par_c, par_zeta, par_kappa, par_n_f, par_n_d, par_epsilon, par_t_n, par_t_A)
% --------------------------------------------------------------------------
function hessp = hessiansp(t, kmrgd, par_C0, par_c, par_zeta, par_kappa, par_n_f, par_n_d, par_epsilon, par_t_n, par_t_A)
%---------------------------------------------------------------------------
function tens3  = der3(t, kmrgd, par_C0, par_c, par_zeta, par_kappa, par_n_f, par_n_d, par_epsilon, par_t_n, par_t_A)
%---------------------------------------------------------------------------
function tens4  = der4(t, kmrgd, par_C0, par_c, par_zeta, par_kappa, par_n_f, par_n_d, par_epsilon, par_t_n, par_t_A)
%---------------------------------------------------------------------------
function tens5  = der5(t, kmrgd, par_C0, par_c, par_zeta, par_kappa, par_n_f, par_n_d, par_epsilon, par_t_n, par_t_A)