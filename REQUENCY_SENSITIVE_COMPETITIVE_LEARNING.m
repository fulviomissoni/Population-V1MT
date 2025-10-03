%% FREQUENCY-SENSITIVE COMPETITIVE LEARNING
% This ensures ALL neurons learn by forcing balanced competition

clear; close all; clc;
addpath FUNCTIONS

%% SETUP (same as before)
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
eta_init = 0.5;
eta_final = 0.05;
max_epochs = 300;  % More epochs for balanced learning

%% RANDOM INITIALIZATION
fprintf('Initializing...\n');
W = randn(n_neurons, n_components) * 0.5 + 0.5;
for i = 1:n_neurons
    W(i,:) = W(i,:) / norm(W(i,:));
end
W_initial = W;

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

fprintf('Mean stimulus correlation: %.3f\n', ...
    mean(corrcoef(activities_all'), 'all'));

%% FREQUENCY-SENSITIVE LEARNING
fprintf('\nStarting frequency-sensitive learning...\n');

winner_counts = zeros(n_neurons, 1);
weight_changes = zeros(max_epochs, 1);

% Track which neuron should get priority
neuron_priorities = ones(n_neurons, 1);

for epoch = 1:max_epochs
    eta = eta_init + (eta_final - eta_init) * (epoch / max_epochs);
    perm = randperm(n_train);
    epoch_change = 0;
    
    for i = 1:n_train
        x = activities_all(:, perm(i));
        
        % Compute distances
        distances = sum((W - x').^2, 2);
        
        % BIAS toward under-utilized neurons
        % Calculate how much each neuron has won relative to fair share
        total_presentations = (epoch - 1) * n_train + i;
        fair_share = total_presentations / n_neurons;
        
        % Neurons below fair share get distance reduction (easier to win)
        for n = 1:n_neurons
            deficit = fair_share - winner_counts(n);
            if deficit > 0
                % Reduce distance proportional to deficit
                bias_factor = 1 - min(0.5, deficit / fair_share);
                distances(n) = distances(n) * bias_factor;
            end
        end
        
        % Find winner with bias applied
        [~, winner] = min(distances);
        winner_counts(winner) = winner_counts(winner) + 1;
        
        % Update winner
        W_old = W(winner, :);
        W(winner, :) = W(winner, :) + eta * (x' - W(winner, :));
        W(winner, :) = W(winner, :) / norm(W(winner, :));
        
        epoch_change = epoch_change + norm(W(winner, :) - W_old);
    end
    
    weight_changes(epoch) = epoch_change / n_train;
    
    if mod(epoch, 25) == 0
        active = sum(winner_counts > 0);
        min_wins = min(winner_counts);
        max_wins = max(winner_counts);
        fprintf('Epoch %d: Active=%d/%d, Wins=[%d-%d], ΔW=%.4f\n', ...
            epoch, active, n_neurons, min_wins, max_wins, weight_changes(epoch));
    end
end

%% POST-TRAINING REFINEMENT FOR DEAD NEURONS
% Force remaining dead neurons to learn from underrepresented stimuli
dead_neurons = find(winner_counts == 0);
if ~isempty(dead_neurons)
    fprintf('\nRefining %d dead neurons...\n', length(dead_neurons));
    
    % For each dead neuron, assign it to most different stimulus
    for d = 1:length(dead_neurons)
        neuron_idx = dead_neurons(d);
        
        % Find which stimuli are least represented by current winners
        stim_coverage = zeros(n_train, 1);
        for s = 1:n_train
            x = activities_all(:, s);
            distances = sum((W - x').^2, 2);
            stim_coverage(s) = min(distances);  % Distance to nearest neuron
        end
        
        % Train this neuron on the worst-covered stimulus
        [~, worst_stim] = max(stim_coverage);
        x_train = activities_all(:, worst_stim);
        
        % Multiple training iterations
        for iter = 1:50
            W(neuron_idx, :) = W(neuron_idx, :) + 0.3 * (x_train' - W(neuron_idx, :));
            W(neuron_idx, :) = W(neuron_idx, :) / norm(W(neuron_idx, :));
        end
        
        fprintf('  Neuron %d trained on stimulus %d\n', neuron_idx, worst_stim);
    end
end

%% ANALYSIS
fprintf('\n=== FINAL RESULTS ===\n');
active_neurons = sum(winner_counts > 0);
fprintf('Active neurons: %d/%d (%.1f%%)\n', active_neurons, n_neurons, ...
    100*active_neurons/n_neurons);
fprintf('Winner distribution: min=%d, max=%d, std=%.1f\n', ...
    min(winner_counts), max(winner_counts), std(winner_counts));

% Compute selectivity for each neuron
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
test_vels = [-4, 0, 4];
test_dirs = [0, pi/2];

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
figure('Position', [100, 100, 1400, 900]);

subplot(2,5,1);
plot(weight_changes, 'LineWidth', 2);
xlabel('Epoch'); ylabel('Weight Change');
title('Learning'); grid on;

subplot(2,5,2);
bar(winner_counts);
xlabel('Neuron'); ylabel('Wins');
title('Winner Distribution'); grid on;

subplot(2,5,3);
histogram(selectivity, 15);
xlabel('Selectivity'); ylabel('Count');
title('Selectivity Distribution'); grid on;

subplot(2,5,4);
imagesc(W_initial);
colorbar; title('Initial');

subplot(2,5,5);
imagesc(W);
colorbar; title('Final');

% Tuning curves
for i = 1:min(20, n_neurons)
    subplot(5, 4, 5+i-1);
    w_2d = reshape(W(i, :), [n_orient, n_vel]);
    imagesc(param.pref_vel, 1:n_orient, w_2d);
    colorbar;
    title(sprintf('N%d w=%d s=%.2f', i, winner_counts(i), selectivity(i)));
    if i > 16, xlabel('Vel'); end
    if mod(i-1, 4) == 0, ylabel('Or'); end
end

%% SAVE
save_filename = sprintf('freq_sensitive_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
save(save_filename, 'W', 'winner_counts', 'selectivity', 'param');
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