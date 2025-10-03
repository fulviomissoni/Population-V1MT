%% REVISED HEBBIAN LEARNING - FIXING ZERO-VELOCITY BIAS
% Key changes: Better stimulus generation, stronger diversity enforcement

clear; close all; clc;
addpath FUNCTIONS

%% PARAMETERS (unchanged)
filter_file = 'FILTERS/New/Gt43B0.0210f0.063.mat';
k0 = 0.063;
samples = 127;
n_orient = 8;
filter_sample = 43;
ft_choice = 'gabor';
v = linspace(-1,1,11)*round(1/2*1/k0);
sigma_pool = 3;
alpha1 = 1e-17;
alpha2 = 1;

param.spat_freq = k0;
param.n_orient = n_orient;
param.pref_vel = v;
param.temp_filt = ft_choice;
param.spatial_filt = filter_file;
param.samples = samples;
param.filter_sample = filter_sample;
param.norm_param = [alpha1, alpha2];
param.sigma_pool = sigma_pool;
param.num_or_ch_pooled = 8;

sigma_orient = 10;
x = linspace(1,n_orient,100);
x_8 = round(linspace(1,100,n_orient));
w_o = exp(-(x(x_8)-(1:n_orient)').^2./(2.*sigma_orient.^2));
param.orient_weighting = w_o./sum(w_o(:));

%% REVISED LEARNING PARAMETERS
learning_params = struct();
learning_params.eta = 0.08;              % Higher learning rate
learning_params.max_epochs = 300;
learning_params.decay_rate = 0.995;
learning_params.k_winners = 3;           % FEWER winners = stronger competition
learning_params.leaky_eta = 0.0005;
learning_params.diversity_penalty = 0.5; % STRONGER penalty for repeated winners

n_pattern_neurons = 25;
n_component_neurons = length(param.pref_vel) * param.n_orient;

%% CRITICAL FIX 1: STRONGER INITIAL DIVERSITY
fprintf('Initializing with STRONG diversity...\n');

W_MT = zeros(n_pattern_neurons, n_component_neurons);
n_vel = length(param.pref_vel);

% Assign each neuron VERY different preferences
vel_assignments = round(linspace(1, n_vel, n_pattern_neurons));
orient_assignments = round(linspace(1, n_orient, n_pattern_neurons));

for i = 1:n_pattern_neurons
    pref_vel_idx = vel_assignments(i);
    pref_orient_idx = orient_assignments(i);
    
    % Create SHARP tuning (narrower than before)
    for vel_idx = 1:n_vel
        for or_idx = 1:n_orient
            comp_idx = (vel_idx-1)*n_orient + or_idx;
            
            vel_dist = abs(vel_idx - pref_vel_idx);
            or_dist = min(abs(or_idx - pref_orient_idx), n_orient - abs(or_idx - pref_orient_idx));
            
            % MUCH sharper tuning (sigma reduced from 4 to 2)
            W_MT(i, comp_idx) = exp(-2 * (vel_dist^2 + or_dist^2));
        end
    end
    
    % Add strong random component to break symmetry
    W_MT(i, :) = W_MT(i, :) + 0.5 * rand(1, n_component_neurons);
    W_MT(i, :) = max(0, W_MT(i, :));
    W_MT(i, :) = W_MT(i, :) / (norm(W_MT(i, :)) + eps);
end

W_MT_initial = W_MT;

%% CRITICAL FIX 2: BETTER STIMULUS GENERATION
fprintf('Generating DIVERSE stimuli (avoiding zero velocity)...\n');

% EXCLUDE zero velocity, use asymmetric range
vel_negative = linspace(-6, -1, 4);  % Negative velocities
vel_positive = linspace(1, 6, 4);    % Positive velocities
vel_range = [vel_negative, vel_positive];  % NO ZERO!

dir_range = linspace(0, 2*pi-2*pi/16, 16);  % Full circle
coh_range = [0.6, 0.8];  % Higher coherence = clearer motion

[V, D, C] = ndgrid(vel_range, dir_range, coh_range);
V = V(:); D = D(:); C = C(:);
n_train_stimuli = length(V);

fprintf('Total training stimuli: %d (no zero velocity)\n', n_train_stimuli);

% Generate stimuli with validation
component_activities_cache = zeros(n_component_neurons, n_train_stimuli);
valid_stimuli = true(n_train_stimuli, 1);

for stim_idx = 1:n_train_stimuli
    stim = create_moving_RDS(V(stim_idx), D(stim_idx), C(stim_idx), param);
    C1 = popFlowV1MT(stim.image, param);
    
    activities = extract_component_activities(C1, param);
    
    % VALIDATION: Check if stimulus produced meaningful response
    if std(activities) < 1e-4 || max(activities) < 1e-6
        warning('Stimulus %d produced weak response, marking invalid', stim_idx);
        valid_stimuli(stim_idx) = false;
    end
    
    component_activities_cache(:, stim_idx) = activities;
    
    if mod(stim_idx, 40) == 0
        fprintf('  Generated %d/%d stimuli\n', stim_idx, n_train_stimuli);
    end
end

% Remove invalid stimuli
V = V(valid_stimuli);
D = D(valid_stimuli);
C = C(valid_stimuli);
component_activities_cache = component_activities_cache(:, valid_stimuli);
n_train_stimuli = sum(valid_stimuli);

fprintf('Valid stimuli: %d\n', n_train_stimuli);

% DIAGNOSTIC: Check stimulus diversity
fprintf('\n=== STIMULUS DIVERSITY CHECK ===\n');
activity_corr = corrcoef(component_activities_cache');
mean_corr = mean(activity_corr(triu(true(size(activity_corr)), 1)));
fprintf('Mean correlation between stimuli: %.3f (should be < 0.7)\n', mean_corr);

if mean_corr > 0.8
    warning('Stimuli are too similar! This will prevent learning.');
end

%% CRITICAL FIX 3: STRONGER COMPETITION
fprintf('\nStarting learning with STRONG competition...\n');

eta = learning_params.eta;
k_winners = learning_params.k_winners;

weight_changes = zeros(learning_params.max_epochs, 1);
winner_diversity = zeros(learning_params.max_epochs, 1);
winner_entropy = zeros(learning_params.max_epochs, 1);

% Win history for diversity enforcement
win_history = zeros(n_pattern_neurons, 1);
consecutive_wins = zeros(n_pattern_neurons, 1);

for epoch = 1:learning_params.max_epochs
    epoch_weight_change = 0;
    winners_count = zeros(n_pattern_neurons, 1);
    
    perm = randperm(n_train_stimuli);
    
    for i = 1:n_train_stimuli
        stim_idx = perm(i);
        component_activities = component_activities_cache(:, stim_idx);
        
        % Compute responses
        pattern_responses = W_MT * component_activities;
        
        % STRONG diversity enforcement
        win_penalty = (win_history / (epoch + 1)).^learning_params.diversity_penalty;
        adjusted_responses = pattern_responses .* (1 - 0.7 * win_penalty);
        
        % Find top-k winners
        [sorted_resp, sorted_idx] = sort(adjusted_responses, 'descend');
        top_k_winners = sorted_idx(1:k_winners);
        
        winners_count(top_k_winners) = winners_count(top_k_winners) + 1;
        win_history(top_k_winners) = win_history(top_k_winners) + 1;
        
        % Update winners with STRONG Hebbian rule
        for k = 1:k_winners
            winner_idx = top_k_winners(k);
            W_before = W_MT(winner_idx, :);
            
            % Simpler, stronger Hebbian update
            rank_factor = 1.0 / k;
            delta_W = eta * rank_factor * (component_activities' - W_MT(winner_idx, :));
            
            W_MT(winner_idx, :) = W_MT(winner_idx, :) + delta_W;
            epoch_weight_change = epoch_weight_change + norm(W_MT(winner_idx, :) - W_before);
        end
        
        % Minimal leaky learning for non-winners
        non_winners = sorted_idx(k_winners+1:end);
        if learning_params.leaky_eta > 0
            for nw_idx = non_winners'
                if pattern_responses(nw_idx) > 0.01
                    W_MT(nw_idx, :) = W_MT(nw_idx, :) + ...
                        learning_params.leaky_eta * (component_activities' - W_MT(nw_idx, :));
                end
            end
        end
        
        % Normalize weights
        for neuron = 1:n_pattern_neurons
            W_MT(neuron, :) = W_MT(neuron, :) / (norm(W_MT(neuron, :)) + eps);
        end
    end
    
    eta = eta * learning_params.decay_rate;
    
    % Compute metrics
    winner_prob = winners_count / sum(winners_count);
    winner_prob = winner_prob(winner_prob > 0);
    
    weight_changes(epoch) = epoch_weight_change / (n_train_stimuli * k_winners);
    winner_diversity(epoch) = sum(winners_count > 0);
    winner_entropy(epoch) = -sum(winner_prob .* log(winner_prob + eps));
    
    if mod(epoch, 20) == 0
        fprintf('Epoch %3d: ΔW=%.4f, Active=%d/%d, Entropy=%.2f/%.2f\n', ...
            epoch, weight_changes(epoch), winner_diversity(epoch), ...
            n_pattern_neurons, winner_entropy(epoch), log(n_pattern_neurons));
    end
    
    if epoch > 50 && std(weight_changes(epoch-9:epoch)) < 1e-5
        fprintf('Converged at epoch %d\n', epoch);
        break;
    end
end

final_epoch = epoch;

%% ANALYSIS
fprintf('\n=== RESULTS ===\n');
fprintf('Final active neurons: %d/%d (%.1f%%)\n', ...
    winner_diversity(final_epoch), n_pattern_neurons, ...
    100*winner_diversity(final_epoch)/n_pattern_neurons);

% Test on diverse velocities (excluding zero)
test_vels = [-5, -3, 3, 5];
test_dirs = [0, pi/4, pi/2, 3*pi/4];

fprintf('\nTesting selectivity:\n');
for v_idx = 1:length(test_vels)
    for d_idx = 1:length(test_dirs)
        test_stim = create_moving_RDS(test_vels(v_idx), test_dirs(d_idx), 0.7, param);
        C1_test = popFlowV1MT(test_stim.image, param);
        act_test = extract_component_activities(C1_test, param);
        resp_test = W_MT * act_test;
        
        [~, best] = maxk(resp_test, 3);
        fprintf('  V=%+.1f, θ=%.2f → [%s]\n', ...
            test_vels(v_idx), test_dirs(d_idx), num2str(best'));
    end
end

%% VISUALIZATION
figure('Position', [100, 100, 1400, 900]);

subplot(2,4,1);
plot(1:final_epoch, weight_changes(1:final_epoch), 'LineWidth', 2);
xlabel('Epoch'); ylabel('Weight Change');
title('Learning Progress'); grid on;

subplot(2,4,2);
plot(1:final_epoch, winner_diversity(1:final_epoch), 'LineWidth', 2);
xlabel('Epoch'); ylabel('Active Neurons');
title('Diversity'); grid on;
ylim([0, n_pattern_neurons]);

subplot(2,4,3);
plot(1:final_epoch, winner_entropy(1:final_epoch), 'LineWidth', 2);
hold on;
plot([1, final_epoch], [log(n_pattern_neurons), log(n_pattern_neurons)], 'r--');
xlabel('Epoch'); ylabel('Entropy');
title('Competition Balance'); grid on;

subplot(2,4,4);
bar(winners_count);
xlabel('Neuron'); ylabel('Wins');
title('Winner Distribution'); grid on;

subplot(2,4,5);
imagesc(W_MT_initial);
colorbar; title('Initial Weights');

subplot(2,4,6);
imagesc(W_MT);
colorbar; title('Final Weights');

subplot(2,4,7);
imagesc(W_MT - W_MT_initial);
colorbar; title('Weight Change');

subplot(2,4,8);
imagesc(activity_corr);
colorbar; caxis([0 1]);
title('Stimulus Correlation Matrix');

% Individual neuron tuning
figure('Position', [100, 100, 1200, 900]);
for i = 1:min(25, n_pattern_neurons)
    subplot(5, 5, i);
    w_2d = reshape(W_MT(i, :), [n_orient, n_vel]);
    imagesc(param.pref_vel, 1:n_orient, w_2d);
    colorbar;
    title(sprintf('N%d (w=%d)', i, winners_count(i)));
    xlabel('Velocity'); ylabel('Orientation');
end

%% SAVE
save_filename = sprintf('hebbian_fixed_zerobias_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
save(save_filename, 'W_MT', 'W_MT_initial', 'learning_params', 'param');
fprintf('\nSaved: %s\n', save_filename);

%% DO NOT MODIFY THIS FUNCTION
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