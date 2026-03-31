
% set scriptPath current directory
directory = 'C:\Users\lenovo\Documents\Masterarbeit\Programming_final_folder';
cd(directory);
% add MatCont files
matContPath = fullfile('C:\Users\lenovo\Documents\MATLAB\MatCont7p4');
addpath(genpath(matContPath));

% Create folder to save outputs
outputFolder = [directory '\output\score'];
savePath = fullfile(outputFolder);
%% Initial conditions
% Spherical condition and Arc length derivative weighed by max(A_c)/max(a_M)
h = @(x) x - x.^2;
dh = @(x) sqrt(1 + (2 - 4*x).^2);

% Define the x range
x_values = linspace(0, 1, 10000);  % Fine grid for integration

% Compute arc length at each x' int(dh, 0, x')
arc_length = zeros(size(x_values));
for i = 2:length(x_values)
    arc_length(i) = integral(dh, 0, x_values(i));
end

% Convert to arc length values
aM_low  = 0.50;
aM_high = 0.95;

s_low  = interp1(x_values, arc_length, aM_low);
s_high = interp1(x_values, arc_length, aM_high);

% Create equal arc-length spacing only on this interval
spacings = 21;
s_values = linspace(s_low, s_high, spacings);

% Interpolate back to get a_M values
a_M_vals = interp1(arc_length, x_values, s_values, 'linear');
a_c_vals = h(a_M_vals);

%% Calculate number of trajectories and score of initial states
% Time scales, rescaled by t_cont, since topological transition is in 
% direction of t_cont and changing t_cont would alter running time till
% transition significantly)
running_time = 3;

% Calculate number of trajectories considered and score expected from (relatively large zeta)
traj_number = 0;
abs_score = 0;
% Loop over the grid of initial points
for i = 1:length(a_M_vals)
    a_M_Init = a_M_vals(i);
    a_c_Init = a_c_vals(i);

    % Exclude small vesicles
    traj_number = traj_number+1;
    abs_score = abs_score + a_M_Init^2*(1-a_M_Init); % a_M * n_fI/n_f * n_dII/n_d at uniform distribution
end
disp(abs_score);
disp(traj_number);

%% Calculate score array

% Parameter ranges
zeta_list = linspace(0, 20, 51);
t_n_list = linspace(1/4, 1.5, 51);

% Default parameters
asize = 1/200;
kappa = 20*4*pi*asize;
epsilon = -1;
n = 0.5;
t_a = 1/5;
C0 = 0.4;
c = 0;

result = [];
    
for t_n = t_n_list
    for zeta = zeta_list
        [score_array, error, slow_mode_count] = calculate_score(a_M_vals, a_c_vals, n/2, n/2, epsilon, kappa, zeta, C0, c, t_n, t_a, running_time);
        
        if size(score_array, 1) > 0
            score_sum = sum(score_array(:,1));
        else
            score_sum = 0;
        end
        
        % Store results including param_value
        result = [result; [t_n, zeta, score_sum/abs_score, error/length(a_M_vals), slow_mode_count/length(a_M_vals)]];
    end
end
% Define filename with param_name and param_value
filename = sprintf('Score_t_n.csv');
csvwrite(fullfile(savePath, filename), result);

%% zeta and c
zeta_list = linspace(0, 8, 51);
c_list = linspace(0, 10, 51); %c of order of sqrt(A/a)
n = 0.1;
t_n = 1;

result = [];
    
for c = c_list
    for zeta = zeta_list
        [score_array, error, slow_mode_count] = calculate_score(a_M_vals, a_c_vals, n/2, n/2, epsilon, kappa, zeta, C0, c, t_n, t_a, running_time);
        
        if size(score_array, 1) > 0
            score_sum = sum(score_array(:,1));
        else
            score_sum = 0;
        end
        
        % Store results including param_value
        result = [result; [c, zeta, score_sum/abs_score, error/length(a_M_vals), slow_mode_count/length(a_M_vals)]];
    end
end
% Define filename with param_name and param_value
filename = sprintf('Score_c.csv');
csvwrite(fullfile(savePath, filename), result);

%% n and c

% Parameter ranges
n_list = linspace(0.025, 0.5, 51);
c_list = linspace(0, 15, 51); %c of order of sqrt(A/a)

% Default parameters
asize = 1/200;
kappa = 20*4*pi*asize;
zeta = 1.5;
epsilon = -1;
t_n = 1;
t_a = 1;
c = 0;
C0 = 0.4;


result = [];

for n = n_list
    for c = c_list
        [score_array, error, slow_mode_count] = calculate_score(a_M_vals, a_c_vals, n/2, n/2, epsilon, kappa, zeta, C0, c, t_n, t_a, running_time);
        
        if size(score_array, 1) > 0
            score_sum = sum(score_array(:,1));
        else
            score_sum = 0;
        end
        
        % Store results including param_value
        result = [result; [n, c, score_sum/abs_score, error/length(a_M_vals), slow_mode_count/length(a_M_vals)]];
    end
end
% Define filename with param_name and param_value
filename = sprintf('Score_n_c.csv');
csvwrite(fullfile(savePath, filename), result);

%% Calculation of the score
function [score_array, error, slow_mode_count] = calculate_score(a_M_vals, a_c_vals, n_f, n_d, epsilon, kappa, zeta, C0, c, t_n, t_A, running_time)
    
    % System for ODEs
    MySystem = Dynamic_system;

    % Initalize score  and slow mode and error count
    score_array = [];
    error = 0;
    slow_mode_count = 0;

    % Loop over the grid of initial points
    for i = 1:length(a_M_vals)
        a_M_Init = a_M_vals(i);
        a_c_Init = a_c_vals(i);


        % Set resolution of  ODE and initialize timer for slow modes
        t_start = tic;
        resolution = 1e-10;

        OPTIONS = odeset('Events', @(t, y) eventFunction(t, y, t_start, n_f, n_d), 'RelTol', resolution, 'AbsTol', resolution);

        % Integrate in time, stop if integration takes too long (slow
        % mode)
        [t, X] = ode45(@(t, y) MySystem{2}(t, y, C0, c, zeta, kappa, n_f, n_d, epsilon, t_n, t_A), ...
            [0 running_time], [a_M_Init, a_c_Init, n_f*a_M_Init, n_d*a_M_Init], OPTIONS);

        final_point = X(end, :);

        [vesicleFormation, invalid, slow_mode] = final_state(t_start, final_point(1), final_point(2), final_point(3), final_point(4), n_f, n_d);
        
        if vesicleFormation
            score_val = final_point(1) * final_point(3)/n_f * (n_d-final_point(4))/n_d;
            score_array = [score_array; score_val];
        elseif invalid
            if real(1-final_point(1)-final_point(2)) < 0.0001 % close to boundary 
                slow_mode_count = slow_mode_count+1;
            else
                disp(X(end-5:end, :));
                disp([t_n, c, zeta, n_f]);
                error=error+1;
                if real((n_f+n_d)-(final_point(3)+final_point(4)))>real(1-final_point(1))
                    disp("Too many particles in the vesicle");
                end
            end
        elseif slow_mode
            slow_mode_count = slow_mode_count+1;
        end
    end
end

function [value, isterminal, direction] = eventFunction(t, y, t_start, n_f, n_d)
    a_M = y(1);
    a_c = y(2);
    n_fI = y(3);
    n_dI = y(4);

    % Check for topological boundary reached
    [vesicleFormation, invalid, slow_mode] = final_state(t_start, a_M, a_c, n_fI, n_dI, n_f, n_d);

    % Stop integration if any condition met
    value = ~(vesicleFormation || invalid || slow_mode);  % 0 if any condition met
    isterminal = 1;
    direction = 0;
end

function [vesicleFormation, invalid, slow_mode] = final_state(t_start, a_M, a_c, n_fI, n_dI, n_f, n_d)

    % Initialize output
    vesicleFormation = 0;
    invalid = 0;
    slow_mode = 0;

    % Check for vesicle formation condition
    if a_c < 0.0001 && all(abs([a_M, a_c, n_fI/n_f, n_dI/n_d]) < 1) && all(abs([a_M, a_c, n_fI, n_dI]) > 0)
        vesicleFormation = 1;

    % Check invalid conditions (imaginary values, invalid geometries,
    % maximum densities)
    elseif any(imag([a_M, a_c, n_fI, n_dI]) ~= 0) ...
        || min(a_M, 1-a_M)-a_c < 0  ...
        || min([n_fI, n_dI, 1-n_fI/a_M, 1-n_dI/a_M, ...
        1-(n_f - n_fI)/(1-a_M), 1-(n_d - n_dI)/(1-a_M)]) < 0

        invalid = 1;

    % Integration time takes too long or vesicle too small (protein size 1/200 and embedded in membrane)
    elseif toc(t_start) > 1 || real(1-a_M)<1/100
        slow_mode = 1;
    end
end