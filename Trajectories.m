% set scriptPath current directory
directory = 'C:\Users\lenovo\Documents\Masterarbeit\Programming_final_folder';
cd(directory);
% add MatCont files
matContPath = fullfile('C:\Users\lenovo\Documents\MATLAB\MatCont7p4');
addpath(genpath(matContPath));

% Create folder to save outputs
outputFolder = [directory '\output\trajectories_flow'];
savePath = fullfile(outputFolder);

%% Trajectories

% Load systems
MySystem2D = Dynamic_system_2D;
MySystem4D = Dynamic_system;
% Set standard options for ODE continuation
OPTIONS=[];

% Parameters
kappa = 4*pi*20;
C0 = 0.4;
zeta = 7/3*4/sqrt(pi)*kappa*(1-C0);

global epsilon

epsilon=-1;

% Initial conditions
a_M_Init = 0.51;
a_c_Init = a_M_Init- a_M_Init^2;

%% 2D
[t, X] = ode45(MySystem2D{2}, [0 100], [a_M_Init, a_c_Init], OPTIONS, C0, zeta, kappa, 1/5);

% Stop when a_c reaches zero
stop_idx = find(real(X(:,2)) <= 0, 1);
if ~isempty(stop_idx)
    t = t(1:stop_idx-1);
    X = X(1:stop_idx-1,:);
end

% Combine results into a table
result = [t, real(X)];

% Save to CSV
fileName = sprintf('2D_main.csv');
csvwrite(fullfile(savePath, fileName), result);

% Annex
zeta = 5/4*4/sqrt(pi)*kappa*(1-C0);

[t, X] = ode45(MySystem2D{2}, [0 100], [a_M_Init, a_c_Init], OPTIONS, C0, zeta, kappa, 1/3);

% Stop when a_c reaches zero
stop_idx = find(real(X(:,2)) <= 0, 1);
if ~isempty(stop_idx)
    t = t(1:stop_idx-1);
    X = X(1:stop_idx-1,:);
end

% Combine results into a table
result = [t, real(X)];

% Save to CSV
fileName = sprintf('2D_appendix.csv');
csvwrite(fullfile(savePath, fileName), result);

%% 4D
particle_frac = 1/200; %a/A
kappa = 20*4*pi*particle_frac;
zeta = 7/3*4/sqrt(pi)*kappa*(1-C0);

% 4D parameter sets
param_sets = [
    struct('n_f', 0.25, 'n_d', 0.25, 'c', 0);
    struct('n_f', 1/20, 'n_d', 1/20, 'c', 0);
    struct('n_f', 1/20, 'n_d', 1/20, 'c', 0.45)
];

for ps = param_sets'
    [n_fM_Init, n_dM_Init] = equi(a_M_Init, ps.n_f, ps.n_d);

    [t, X] = ode45(MySystem4D{2}, [0 1000], [a_M_Init, a_c_Init, n_fM_Init, n_dM_Init], OPTIONS, C0, ps.c, zeta, kappa, ps.n_f, ps.n_d, epsilon, 1/1000, 1/5);

    % Stop when a_c = 0
    stop_idx = find(real(X(:,2)) <= 0, 1);
    if ~isempty(stop_idx)
        t = t(1:stop_idx-1);
        X = X(1:stop_idx-1,:);
    end

    % Combine results into a table
    result = [t, real(X)];
    
    % Save to CSV
    fileName = sprintf('4D_c%.2f_nf%.2f_main.csv', ps.c, ps.n_f);;
    csvwrite(fullfile(savePath, fileName), result);
end

% Appendix

zeta = 4/3*4/sqrt(pi)*kappa*(1-C0);

for ps = param_sets'
    [n_fM_Init, n_dM_Init] = equi(a_M_Init, ps.n_f, ps.n_d);

    [t, X] = ode45(MySystem4D{2}, [0 1000], [a_M_Init, a_c_Init, n_fM_Init, n_dM_Init], OPTIONS, C0, ps.c, zeta, kappa, ps.n_f, ps.n_d, epsilon, 1/1000, 1/3);

    % Stop when a_c = 0
    stop_idx = find(real(X(:,2)) <= 0, 1);
    if ~isempty(stop_idx)
        t = t(1:stop_idx-1);
        X = X(1:stop_idx-1,:);
    end

    % Combine results into a table
    result = [t, real(X)];
    
    % Save to CSV
    fileName = sprintf('4D_c%.2f_nf%.2f_appendix.csv', ps.c, ps.n_f);;
    csvwrite(fullfile(savePath, fileName), result);
end

%% Find equilibrium for particle numbers
% Equilibrium conditions for c = 0

function [n_fM, n_dM] = equi(a_M, n_f, n_d)
    global	epsilon

    gamma = 1 / (1 + exp(-epsilon));

    n_fM = ((a_M + n_f) * (1 - gamma) + gamma - ...
        sqrt(-4 * a_M * n_f * (1 - gamma) + ((a_M + n_f) * (1 - gamma) + gamma).^2)) / (2 * (1 - gamma));

    n_dM = (-a_M * n_d + (a_M * n_d) / (2 * (1 - gamma)) ...
        + (n_f * n_d) / (2 * (1 - gamma)) + (n_d * gamma) / (2 * (1 - gamma)) ...
        - (a_M * n_d * gamma) / (2 * (1 - gamma)) - (n_f * n_d * gamma) / (2 * (1 - gamma)) ...
        - (n_d * sqrt(-4 * a_M * n_f * (1 - gamma) ...
        + (a_M + n_f + gamma - a_M * gamma - n_f * gamma).^2)) / (2 * (1 - gamma))) / (n_f - 1);
end