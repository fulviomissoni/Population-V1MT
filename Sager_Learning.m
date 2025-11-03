%% COMPETITIVE LEARNING WITH HOMEOSTATIC PLASTICITY
% Biologically plausible: neurons that lose consistently become more excitable

clear; close all; clc;
addpath FUNCTIONS

%% SETUP
filter_file = 'FILTERS/New/Gt43B0.0210f0.063.mat';
k0 = 0.063; samples = 127; n_orient = 8; filter_sample = 43;
ft_choice = 'gabor';
v = linspace(-1,1,11)*round(1/2*1/k0);

param.spat_freq = k0; param.n_orient = n_orient; param.pref_vel = v;
param.temp_filt = ft_choice; param.spatial_filt = filter_file;
param.samples = samples; param.filter_sample = filter_sample;
param.norm_param = [1e-17, 1]; param.sigma_pool = 3; param.num_or_ch_pooled = 8;

sigma_orient = 10; x = linspace(1,n_orient,100); x_8 = round(linspace(1,100,n_orient));
w_o = exp(-(x(x_8)-(1:n_orient)').^2./(2.*sigma_orient.^2));
param.orient_weighting = w_o./sum(w_o(:));

%% PARAMETERS
n_neurons = 25;
n_components = length(param.pref_vel) * param.n_orient;

% Learning parameters
eta_init = 0.5;
eta_final = 0.05;
max_epochs = 300;

% Homeostatic parameters
tau_homeo = 50;           % Time constant for homeostasis (in epochs)
target_rate = 1.0;        % Target average activity level
homeo_strength = 0.02;    % How strongly homeostasis affects excitability

%% RANDOM INITIALIZATION
fprintf('Random initialization...\n');
W = randn(n_neurons, n_components) * 0.5 + 0.5;
for i = 1:n_neurons
    W(i,:) = W(i,:) / norm(W(i,:));
end
W_initial = W;

% Initialize homeostatic variables
excitability = ones(n_neurons, 1);  % Multiplicative gain for each neuron
avg_activity = zeros(n_neurons, 1); % Running average of activity

%% GENERATE TRAINING DATA
fprintf('Generating training data...\n');
vel_vals = param.pref_vel(2:end-1);
dir_vals = linspace(0, 2*pi-2*pi/12, 12);
[VV, DD] = meshgrid(vel_vals, dir_vals);
VV = VV(:); DD = DD(:);
n_train = length(VV);

activities_all = zeros(n_components, n_train);
for i = 1:n_train
    stim = create_moving_RDS(VV(i), DD(i), 0.7, param);
    C1 = popFlowV1MT(stim.image, param);
    activities_all(:, i) = extract_component_activities(C1, param);
    if mod(i, 20) == 0
        fprintf('  %d/%d\n', i, n_train);
    end
end

corr_mat = corrcoef(activities_all');
fprintf('Mean stimulus correlation: %.3f\n', mean(corr_mat(triu(true(size(corr_mat)), 1))));

%% COMPETITIVE LEARNING WITH HOMEOSTASIS
fprintf('\nStarting learning with homeostatic plasticity...\n');

winner_counts = zeros(n_neurons, 1);
weight_changes = zeros(max_epochs, 1);
excitability_history = zeros(max_epochs, n_neurons);

for epoch = 1:max_epochs
    eta = eta_init + (eta_final - eta_init) * (epoch / max_epochs);
    
    perm = randperm(n_train);
    epoch_change = 0;
    epoch_activity = zeros(n_neurons, 1);
    
    for i = 1:n_train
        x = activities_all(:, perm(i));
        
        % Compute distances (similarity to input)
        distances = sum((W - x').^2, 2);
        
        % Apply homeostatic excitability modulation
        % Lower distance = more similar = higher response
        % Multiply by excitability: under-active neurons get boosted
        responses = (1 ./ (distances + 0.01)) .* excitability;
        
        % Winner: highest response (considering both similarity and excitability)
        [~, winner] = max(responses);
        
        winner_counts(winner) = winner_counts(winner) + 1;
        epoch_activity(winner) = epoch_activity(winner) + 1;
        
        % Update winner weights (standard competitive learning)
        W_old = W(winner, :);
        W(winner, :) = W(winner, :) + eta * (x' - W(winner, :));
        W(winner, :) = W(winner, :) / norm(W(winner, :));
        
        epoch_change = epoch_change + norm(W(winner, :) - W_old);
    end
    
    % UPDATE HOMEOSTATIC EXCITABILITY (slow, once per epoch)
    % Neurons with low activity increase excitability
    % Neurons with high activity decrease excitability
    epoch_activity = epoch_activity / n_train;  % Normalize by number of stimuli
    
    % Update running average of activity with exponential filter
    alpha_filter = 1 - exp(-1/tau_homeo);
    avg_activity = (1 - alpha_filter) * avg_activity + alpha_filter * epoch_activity;
    
    % Adjust excitability based on deviation from target
    activity_error = target_rate - avg_activity;
    excitability = excitability + homeo_strength * activity_error;
    
    % Keep excitability within reasonable bounds
    excitability = max(0.5, min(excitability, 3.0));
    
    % Store metrics
    weight_changes(epoch) = epoch_change / n_train;
    excitability_history(epoch, :) = excitability;
    
    if mod(epoch, 25) == 0
        active = sum(winner_counts > 0);
        mean_exc = mean(excitability);
        std_exc = std(excitability);
        fprintf('Epoch %d: η=%.3f, ΔW=%.4f, Active=%d/%d, Exc=%.2f±%.2f\n', ...
            epoch, eta, weight_changes(epoch), active, n_neurons, mean_exc, std_exc);
    end
end

%% ANALYSIS
fprintf('\n=== RESULTS ===\n');
active_neurons = sum(winner_counts > 0);
fprintf('Active neurons: %d/%d (%.1f%%)\n', active_neurons, n_neurons, ...
    100*active_neurons/n_neurons);
fprintf('Final excitability: mean=%.2f, range=[%.2f, %.2f]\n', ...
    mean(excitability), min(excitability), max(excitability));
fprintf('Total weight change: %.3f\n', norm(W - W_initial, 'fro'));

% Compute selectivity
selectivity = zeros(n_neurons, 1);
n_vel = length(param.pref_vel);
for neuron = 1:n_neurons
    w_2d = reshape(W(neuron, :), [n_orient, n_vel]);
    max_val = max(w_2d(:));
    mean_val = mean(w_2d(:));
    selectivity(neuron) = (max_val - mean_val) / (max_val + eps);
end

fprintf('Selectivity: mean=%.3f, std=%.3f\n', mean(selectivity), std(selectivity));
fprintf('Selective neurons (>0.3): %d/%d\n', sum(selectivity > 0.3), n_neurons);

%% TEST
fprintf('\nTesting:\n');
test_vels = [-4, -2, 0, 2, 4];
test_dirs = [0, pi/4, pi/2];

for v_idx = 1:length(test_vels)
    for d_idx = 1:length(test_dirs)
        stim_test = create_moving_RDS(test_vels(v_idx), test_dirs(d_idx), 0.7, param);
        C1_test = popFlowV1MT(stim_test.image, param);
        x_test = extract_component_activities(C1_test, param);
        
        distances = sum((W - x_test').^2, 2);
        [~, sorted] = sort(distances);
        
        fprintf('  V=%+.1f, θ=%.2f → [%d, %d, %d]\n', ...
            test_vels(v_idx), test_dirs(d_idx), sorted(1), sorted(2), sorted(3));
    end
end

%% VISUALIZATION
figure('Position', [100, 100, 1600, 900]);

% Learning and homeostasis dynamics
subplot(2,5,1);
plot(weight_changes, 'LineWidth', 2);
xlabel('Epoch'); ylabel('Weight Change');
title('Learning Progress'); grid on;

subplot(2,5,2);
plot(excitability_history, 'LineWidth', 1.5);
xlabel('Epoch'); ylabel('Excitability');
title('Homeostatic Excitability'); grid on;
legend('Location', 'best');

subplot(2,5,3);
bar(winner_counts);
xlabel('Neuron'); ylabel('Wins');
title(sprintf('Winners (Active: %d/%d)', active_neurons, n_neurons));
grid on;

subplot(2,5,4);
bar(excitability);
xlabel('Neuron'); ylabel('Final Excitability');
title('Homeostatic Gains'); grid on;
yline(1.0, 'r--', 'Baseline');

subplot(2,5,5);
scatter(winner_counts, excitability, 50, 'filled');
xlabel('Win Count'); ylabel('Excitability');
title('Wins vs Excitability'); grid on;

% Individual neuron tuning
for i = 1:min(20, n_neurons)
    subplot(5, 5, 5+i-1);
    w_2d = reshape(W(i, :), [n_orient, n_vel]);
    imagesc(param.pref_vel, 1:n_orient, w_2d);
    colorbar;
    title(sprintf('N%d w=%d e=%.2f', i, winner_counts(i), excitability(i)));
    if i > 16, xlabel('Vel'); end
    if mod(i-1, 4) == 0, ylabel('Or'); end
end

% Additional analysis figure
figure('Position', [100, 100, 1200, 400]);

subplot(1,3,1);
histogram(selectivity, 15);
xlabel('Selectivity'); ylabel('Count');
title('Selectivity Distribution'); grid on;

subplot(1,3,2);
scatter(excitability, selectivity, 50, winner_counts, 'filled');
colorbar; xlabel('Excitability'); ylabel('Selectivity');
title('Excitability vs Selectivity');
grid on;

subplot(1,3,3);
imagesc([W_initial, W]);
colorbar; 
title('Initial (left) vs Final (right) Weights');
xlabel('Neuron'); ylabel('Component');

figure, popresponse_tiled(permute(reshape(W,5,5,8,11),[3,4,1,2])), title('Learned Weights per MT cell')
figure, popresponse_tiled(reshape(activities_all,8,11,9,12)), title('V1 activity')
%% SAVE
save_filename = sprintf('homeostatic_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
save(save_filename, 'W', 'W_initial', 'excitability', 'winner_counts', ...
    'selectivity', 'excitability_history', 'param');
fprintf('\nSaved: %s\n', save_filename);

%% HELPER
function activities = extract_component_activities(C1, param)
    [sy, sx, n_orient, n_vel] = size(C1{1});
    combined = (C1{2});
    activities = squeeze(mean(combined, [1, 2]));
    activities = activities(:);
    activities = max(0, activities);
    activities = activities + 0.001 * max(activities);
    if norm(activities) > 0
        activities = activities / norm(activities);
    end
end