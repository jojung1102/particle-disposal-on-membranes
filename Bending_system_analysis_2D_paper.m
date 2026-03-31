% set scriptPath current directory
directory = 'C:\Users\lenovo\Documents\Masterarbeit\Programming_final_folder';
cd(directory);
% add MatCont files
matContPath = fullfile('C:\Users\lenovo\Documents\MATLAB\MatCont7p4');
addpath(genpath(matContPath));

% Create folder to save outputs
outputFolder = [directory '\output\Bifurcations'];
savePath = fullfile(outputFolder);

% Load system
MySystem = Bending_system_2D;

% Set standard options for ODE continuation
OPTIONS=[];

% Initial conditions
a_M_Init = 0.5;
a_c_Init = 0.4;
kappa = 20*4*pi;
zeta = 0.5;
C0 = 0.4;

% Integrate in time
[t, X_pt] = ode45(MySystem{2}, [0 1000], [a_M_Init, a_c_Init], OPTIONS, C0, zeta, kappa);

% Initialize equilibrium continuation from the last point of the trajectory
% Branching points are found when entropic contribution is dropped; For
% the main branch, the entropic contribution doesn't play a role
% (nullclines far apart)
p = [C0; zeta; kappa];
ap = [2];  % Active parameter in continuation

% Initial state values (vertical vector)
X1 = [X_pt(end, 1); X_pt(end, 2)];
[X0, ~] = init_EP_EP(@Bending_system_2D, X1, p, ap);

% Options for continuation
opt = contset;
opt = contset(opt, 'MaxNumPoints', 5*1e5);
opt = contset(opt, 'InitStepsize', 0.0001);
opt = contset(opt, 'MaxStepsize', 0.005);
opt = contset(opt, 'Singularities', 1);
opt = contset(opt, 'Eigenvalues', 1);
opt = contset(opt, 'Backward', 0);

% Continuation
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
X_combined = [X_combined; sum(F_combined <= 0, 1)]; % stability: 0 unstable, 1 saddle, 2 stable

% Identify special points
special_points = vertcat(S_back, S);
limit_indices = [special_points(strcmp({special_points.label}, 'LP')).index];
branching_indices = [special_points(strcmp({special_points.label}, 'BP')).index];
branching_points = special_points(strcmp({special_points.label}, 'BP'));

% Label special points
X_combined(4, branching_indices) = 3; % 3 for branching points
X_combined(4, limit_indices) = 4;  % 4 for limit points

% Save combined branch
X_combined_positive = real(X_combined(:, X_combined(3, :) >= 0));
csvwrite(fullfile(savePath, 'Main_branch.csv'), X_combined_positive');

branch = 1;

%% Side branches
% Branch switching from BP
for i = 1:length(branching_points)
    p(ap) = X_combined(3, branching_points(i).index);
    
    [X0_side, ~] = init_BP_EP(@Bending_system_2D, X_combined(1:2, branching_points(i).index), p, branching_points(i), 0.00001);

    opt = contset(opt, 'Backward', 0);
    [X_side, ~, S_side, ~, F_side] = cont(@equilibrium, X0_side, [], opt);

    opt = contset(opt, 'Backward', 1);
    [X_side_back, ~, S_side_back, ~, F_side_back] = cont(@equilibrium, X0_side, [], opt);

    % Reverse columns; F sometimes contains one additional entry at the end point
    difference = length(F_side) - length(X_side);
    disp(difference);
    F_side = F_side(:, 1:end - difference);

    difference = length(F_side_back) - length(X_side_back);
    disp(difference);
    X_side_back_reversed = fliplr(X_side_back);
    F_side_back_reversed = fliplr(F_side_back(:, 1:end - difference));

    % Adjust the indices in S_back
    len_X_side_back = size(X_side_back, 2);
    len_X_side = size(X_side, 2);

    % Update the indices in S
    for j = 1:length(S_side_back)
        S_side_back(j).index = len_X_side_back - S_side_back(j).index + 1;
    end

    for j = 1:length(S_side)
        S_side(j).index = len_X_side_back + S_side(j).index;
    end

    % Combine the two parts of the main branch
    X_side_combined = [X_side_back_reversed, X_side];
    F_side_combined = [F_side_back_reversed, F_side];

    % Main branch with stability
    X_side_combined = [X_side_combined; sum(F_side_combined <= 0, 1)]; % stability: 0 unstable, 1 saddle, 2 stable

    % Identify special points
    special_points_side = vertcat(S_side_back, S_side);
    branching_indices_side = [special_points_side(strcmp({special_points_side.label}, 'BP')).index];
    limit_indices_side = [special_points_side(strcmp({special_points_side.label}, 'LP')).index];

    % Label special points
    X_side_combined(4, branching_indices_side) = 3;  % 3 for branching points
    X_side_combined(4, limit_indices_side) = 4;  % 4 for limit points
    
    % Save branch
    X_side_combined_positive = real(X_side_combined(:, X_side_combined(3, :) >= 0));
    csvwrite(fullfile(savePath, sprintf('Side_branch_%d.csv', branch)), X_side_combined_positive');
    
    branch = branch + 1;
end
