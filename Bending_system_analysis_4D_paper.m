% set scriptPath current directory
directory = 'C:\Users\lenovo\Documents\Masterarbeit\Programming_final_folder';
cd(directory);
% add MatCont files
matContPath = fullfile('C:\Users\lenovo\Documents\MATLAB\MatCont7p4');
addpath(genpath(matContPath));

% Create folder to save outputs
inputFolder = [directory '\output\Bifurcations'];
outputFolder = [directory '\output\Bifurcations'];
inputPath = fullfile(inputFolder);
savePath = fullfile(outputFolder);

% Load system
MySystem = Bending_system_4D;
% Set standard options for ODE continuation
OPTIONS=[];

% Initial conditions t_n, t_a, zeta, kappa are rescaled by a/A
global kappa C0 n n_f n_d e
particle_frac = 1/200; %a/A
kappa = 20*4*pi*particle_frac;  % For more realistic values, the branches change their position too much such that they are not found anymore
C0 = 0.4;
n = 0.1;
n_f = n/2;
n_d = n/2;
e = -1;

% Import main branch
mainBranch = table2array(readtable(fullfile(inputPath, 'Main_branch.csv')));

% Find BP indices
splitIndices = find(mainBranch(:,4) == 3);
splitIndices = [splitIndices; size(mainBranch, 1)];

% Split into sections based at BPs and choose a point on this branch
selectedRows = [];
startIdx = 1;
for i = 1:length(splitIndices)

    endIdx = splitIndices(i);
    section = mainBranch(startIdx:endIdx, :);
    
    % Compute middle of range in column 3
    middleValue = (min(section(:,3)) + max(section(:,3))) / 2;
    
    % Find row closest to the middle value in column 3
    [~, closestIdx] = min(abs(section(:,3) - middleValue));
    
    % Store the selected row
    selectedRows = [selectedRows; section(closestIdx, :)];
    
    % Update start index for the next section
    startIdx = endIdx + 1;
end

%% c=0
branch = 1;
for i = 1:size(selectedRows, 1)
    X_combined = branch_continuation(0, selectedRows(i, 1), selectedRows(i, 2), particle_frac*selectedRows(i, 3));
    csvwrite(fullfile(savePath, sprintf('Branch_c_0_%d.csv', branch)), X_combined');
    branch = branch + 1;
end

%% c=0.45
branch = 1;
for i = 1:size(selectedRows, 1)
    X_combined = branch_continuation(0.45, selectedRows(i, 1), selectedRows(i, 2), particle_frac*selectedRows(i, 3));
    csvwrite(fullfile(savePath, sprintf('Branch_c_0.45_%d.csv', branch)), X_combined');
    branch = branch + 1;
end

%% Functions
% Equilibrium conditions for c = 0
function [n_fM, n_dM] = equi(a_M)
    global	e n_f n_d

    gamma = 1 / (1 + exp(-e));

    n_fM = ((a_M + n_f) * (1 - gamma) + gamma - ...
        sqrt(-4 * a_M * n_f * (1 - gamma) + ((a_M + n_f) * (1 - gamma) + gamma).^2)) / (2 * (1 - gamma));

    n_dM = (-a_M * n_d + (a_M * n_d) / (2 * (1 - gamma)) ...
        + (n_f * n_d) / (2 * (1 - gamma)) + (n_d * gamma) / (2 * (1 - gamma)) ...
        - (a_M * n_d * gamma) / (2 * (1 - gamma)) - (n_f * n_d * gamma) / (2 * (1 - gamma)) ...
        - (n_d * sqrt(-4 * a_M * n_f * (1 - gamma) ...
        + (a_M + n_f + gamma - a_M * gamma - n_f * gamma).^2)) / (2 * (1 - gamma))) / (n_f - 1);
end

% Continuation in c, then in zeta
function X_combined = branch_continuation(c, a_M_Init, a_c_Init, zeta)
    global kappa C0 n_f n_d e
    % Initial variables & parameters
    [n_fM_Init, n_dM_Init] = equi(a_M_Init);
    p = [C0; 0; zeta; kappa; n_f; n_d; e]; % Parameters
    ap = [2];  % Active parameter in continuation
    
    % Initial state values (vertical vector)
    X1 = [a_M_Init; a_c_Init; n_fM_Init; n_dM_Init];
    [X0, ~] = init_EP_EP(@Bending_system_4D, X1, p, ap);
    
    % Options for continuation
    opt = contset;
    opt = contset(opt, 'MaxNumPoints', 1e4);
    opt = contset(opt, 'InitStepsize', 0.001);
    opt = contset(opt, 'MaxStepsize', 0.01); % Rather refine option??
    opt = contset(opt, 'Singularities', 1);
    opt = contset(opt, 'Eigenvalues', 1);
    opt = contset(opt, 'Backward', 0);
    
    % Continuations
    [X, ~, ~, ~, ~] = cont(@equilibrium, X0, [], opt);
    
    opt = contset(opt, 'Backward', 1);
    [X_back, ~, ~, ~, ~] = cont(@equilibrium, X0, [], opt);

    X_back_reversed = fliplr(X_back);

    % Combine the two parts of the main branch
    X_combined = [X_back_reversed, X];

    % Select initial variables and parameters for continuation in zeta
    [~, idx] = min(abs(X_combined(5,:) - c));
    a_M_Init = X_combined(1, idx);
    a_c_Init = X_combined(2, idx);
    n_fM_Init= X_combined(3, idx);
    n_dM_Init = X_combined(4, idx);
    p(2) = c;
    ap = [3];  % Active parameter in continuation

    % Initial state values (vertical vector)
    X1 = [a_M_Init; a_c_Init; n_fM_Init; n_dM_Init];
    [X0, ~] = init_EP_EP(@Bending_system_4D, X1, p, ap);
    
    % Continuations
    opt = contset(opt, 'Backward', 0);
    [X, ~, S, ~, F] = cont(@equilibrium, X0, [], opt);
    
    opt = contset(opt, 'Backward', 1);
    [X_back, ~, S_back, ~, F_back] = cont(@equilibrium, X0, [], opt);

    % Reverse columns; F sometimes contains one additional entry at the end point
    difference = length(F) - length(X);
    F = F(:, 1:end - difference);

    difference = length(F_back) - length(X_back);
    X_back_reversed = fliplr(X_back);
    F_back_reversed = fliplr(F_back(:, 1:end - difference));
    
    % Adjust the indices in S_back
    len_X_back = size(X_back, 2);
    len_X = size(X, 2);
    
    % Update the indices in S
    for j = 1:length(S_back)
        S_back(j).index = len_X_back - S_back(j).index + 1;
    end
    
    for j = 1:length(S)
        S(j).index = len_X_back + S(j).index;
    end
    
    % Combine the two parts of the main branch
    X_combined = [X_back_reversed, X];
    F_combined = [F_back_reversed, F];
    
    % Main branch with stability
    X_combined = [X_combined; sum(F_combined <= 0, 1)-2]; % stability: 0 unstable, 1 saddle, 2 stable
    
    % Identify special points
    special_points = vertcat(S_back, S);
    limit_indices = [special_points(strcmp({special_points.label}, 'LP')).index];
    branching_indices = [special_points(strcmp({special_points.label}, 'BP')).index];
    
    % Label special points
    X_combined(6, branching_indices) = 3; % 3 for branching points
    X_combined(6, limit_indices) = 4;  % 4 for limit points
    
    % Save combined branch
    X_combined = real(X_combined(:, real(X_combined(5, :)) >= 0));
end


