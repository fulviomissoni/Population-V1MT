%% FAST FIXED HEBBIAN LEARNING FOR MT MODEL
% Addresses: 1) Single winner problem, 2) Speed issues

clear
close all
clc
addpath FUNCTIONS

%% PARAMETERS (same as before)
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

%% OPTIMIZED LEARNING PARAMETERS
learning_params = struct();
learning_params.eta = 0.02;              % Slower learning rate
learning_params.max_epochs = 300;
learning_params.conv_threshold = 1e-4;
learning_params.decay_rate = 0.995;
learning_params.k_winners = 15;          % More winners
learning_params.inhibition_strength = 0.0; % NO inhibition initially

n_pattern_neurons = 30;  % Fewer neurons = faster + better learning
n_component_neurons = length(param.pref_vel) * param.n_orient;

%% BETTER INITIALIZATION - SIMPLE AND ROBUST
fprintf('Initializing weights with simple random distribution...\n');

% Completely random initialization (most robust)
W_MT = randn(n_pattern_neurons, n_component_neurons);
n_vel = length(param.pref_vel);
% Add some structured diversity
for i = 1:n_pattern_neurons
    % Give each neuron a slight preference
    pref_vel_idx = randi(length(param.pref_vel));
    pref_or_idx = randi(param.n_orient);
    
    % Boost weights around preference
    center_idx = (pref_vel_idx-1)*param.n_orient + pref_or_idx;
    boost_range = max(1, center_idx-3):min(n_component_neurons, center_idx+3);
    W_MT(i, boost_range) = W_MT(i, boost_range) + 0.5 * randn(1, length(boost_range));
    
    % Normalize
    W_MT(i, :) = W_MT(i, :) / (norm(W_MT(i, :)) + eps);
    
    % Ensure no dead neurons by adding small positive bias
    W_MT(i, :) = W_MT(i, :) + 0.1 * rand(1, n_component_neurons);
    W_MT(i, :) = W_MT(i, :) / (norm(W_MT(i, :)) + eps);
end

W_MT_initial = W_MT;
% %% IMPROVED INITIALIZATION - DIVERSE PREFERENCES
% fprintf('Initializing weights with diverse preferences...\n');
% 
% n_vel = length(param.pref_vel);
% vel_prefs = linspace(1, n_vel, n_pattern_neurons);
% dir_prefs = linspace(0, 2*pi, n_pattern_neurons);
% 
% W_MT = zeros(n_pattern_neurons, n_component_neurons);
% 
% for i = 1:n_pattern_neurons
%     % Each neuron has a preferred velocity and direction
%     pref_vel_idx = round(vel_prefs(i));
%     pref_dir_idx = round(dir_prefs(i) / (2*pi) * n_orient);
% 
%     % Create Gaussian tuning around preference
%     for vel_idx = 1:n_vel
%         for or_idx = 1:n_orient
%             comp_idx = (vel_idx-1)*n_orient + or_idx;
% 
%             % Distance in feature space
%             vel_dist = (vel_idx - pref_vel_idx)^2 / n_vel^2;
% 
%             % Circular distance for orientation
%             or_dist_1 = abs(or_idx - pref_dir_idx);
%             or_dist_2 = n_orient - or_dist_1;
%             or_dist = min(or_dist_1, or_dist_2)^2 / n_orient^2;
% 
%             % Gaussian weight
%             W_MT(i, comp_idx) = exp(-5 * (vel_dist + or_dist));
%         end
%     end
% 
%     % Normalize + add noise
%     W_MT(i, :) = W_MT(i, :) / norm(W_MT(i, :));
%     W_MT(i, :) = W_MT(i, :) + randn(1, n_component_neurons) * 0.05;
%     W_MT(i, :) = W_MT(i, :) / norm(W_MT(i, :));
% end

%% FASTER TRAINING DATA GENERATION - FEWER STIMULI
fprintf('Generating training stimuli (reduced set for speed)...\n');

% Reduce number of stimuli for speed
vel_range = linspace(min(param.pref_vel)*0.5, max(param.pref_vel)*0.5, 22);  % 10 -> 6
dir_range = linspace(0, 2*pi, 12);  % 16 -> 12
coh_range = [0.5, 0.7];  % 4 -> 2

[V, D, C] = ndgrid(vel_range, dir_range, coh_range);
V = V(:); D = D(:); C = C(:);
n_train_stimuli = length(V);

training_data = cell(n_train_stimuli, 1);
component_activities_cache = zeros(n_component_neurons, n_train_stimuli);

fprintf('Total training stimuli: %d\n', n_train_stimuli);

for stim_idx = 1:n_train_stimuli
    stim = create_moving_RDS(V(stim_idx), D(stim_idx), C(stim_idx), param);
    C1 = popFlowV1MT(stim.image, param);
    
    % Cache activities to avoid recomputation
    component_activities_cache(:, stim_idx) = extract_component_activities_fixed(C1, param);
    
    if mod(stim_idx, 30) == 0
        fprintf('  Generated %d/%d stimuli\n', stim_idx, n_train_stimuli);
    end
end

%% TOP-K HEBBIAN LEARNING (NOT WINNER-TAKE-ALL)
fprintf('\nStarting TOP-K Hebbian learning...\n');

eta = learning_params.eta;
k_winners = learning_params.k_winners;

weight_changes = zeros(learning_params.max_epochs, 1);
avg_responses = zeros(learning_params.max_epochs, 1);
winner_diversity = zeros(learning_params.max_epochs, 1);
winner_entropy = zeros(learning_params.max_epochs, 1);

W_MT_initial = W_MT;

for epoch = 1:learning_params.max_epochs
    epoch_weight_change = 0;
    epoch_response = 0;
    winners_count = zeros(n_pattern_neurons, 1);
    
    % Shuffle
    perm = randperm(n_train_stimuli);
    
    for i = 1:n_train_stimuli
        stim_idx = perm(i);
        
        % Use cached activities (MUCH FASTER!)
        component_activities = component_activities_cache(:, stim_idx);
        
        % Compute responses
        pattern_responses = W_MT * component_activities;

        % CRITICAL: Apply response threshold - only consider neurons above threshold
        response_threshold = 0.05 * max(pattern_responses);  % 5% of max response
        valid_neurons = pattern_responses > response_threshold;
        
        if sum(valid_neurons) < learning_params.k_winners
            % If too few neurons respond, use top ones anyway
            [~, sorted_idx] = sort(pattern_responses, 'descend');
            top_k_winners = sorted_idx(1:learning_params.k_winners);
        else
            % Only consider neurons above threshold
            valid_responses = pattern_responses(valid_neurons);
            valid_indices = find(valid_neurons);
            [~, sort_idx] = sort(valid_responses, 'descend');
            top_k_winners = valid_indices(sort_idx(1:min(learning_params.k_winners, length(valid_indices))));
        end

        % TOP-K WINNERS instead of single winner
        [sorted_resp, sorted_idx] = sort(pattern_responses, 'descend');
        % top_k_winners = sorted_idx(1:k_winners);

        % Track statistics
        epoch_response = epoch_response + sorted_resp(1);
        winners_count(top_k_winners) = winners_count(top_k_winners) + 1;
        
        % UPDATE EACH TOP-K WINNER (with decreasing strength)
        for k = 1:k_winners
            winner_idx = top_k_winners(k);
            W_before = W_MT(winner_idx, :);
            
            % Learning rate decreases with rank
            rank_factor = 1.0 / k;
            
            % Hebbian update
            delta_W = eta * rank_factor * component_activities' * sorted_resp(k);
            W_MT(winner_idx, :) = W_MT(winner_idx, :) + delta_W;
            
            % Track change
            epoch_weight_change = epoch_weight_change + norm(W_MT(winner_idx, :) - W_before);
        end
        
        % LATERAL INHIBITION - suppress non-winners
        non_winners = sorted_idx(k_winners+1:end);
        for nw_idx = non_winners'
            W_MT(nw_idx, :) = W_MT(nw_idx, :) * (1 - learning_params.inhibition_strength);
        end
    end
    
    % Normalize all weights
    for i = 1:n_pattern_neurons
        W_MT(i, :) = W_MT(i, :) / (norm(W_MT(i, :)) + eps);
    end
    
    % Decay learning rate
    eta = eta * learning_params.decay_rate;
    
    % Compute diversity metrics
    winner_prob = winners_count / sum(winners_count);
    winner_prob = winner_prob(winner_prob > 0);
    
    weight_changes(epoch) = epoch_weight_change / n_train_stimuli;
    avg_responses(epoch) = epoch_response / n_train_stimuli;
    winner_diversity(epoch) = sum(winners_count > 0);
    winner_entropy(epoch) = -sum(winner_prob .* log(winner_prob + eps));
    
    % Print progress
    if mod(epoch, 25) == 0
        fprintf('Epoch %d: Weight Δ=%.5f, Resp=%.3f, Active=%d/%d, Entropy=%.2f\n', ...
            epoch, weight_changes(epoch), avg_responses(epoch), ...
            winner_diversity(epoch), n_pattern_neurons, winner_entropy(epoch));
    end
    
    % Convergence check
    if epoch > 50
        recent_changes = weight_changes(max(1, epoch-10):epoch);
        if std(recent_changes) < learning_params.conv_threshold
            fprintf('Converged at epoch %d\n', epoch);
            break;
        end
    end
end

%% ANALYSIS
fprintf('\n=== RESULTS ===\n');
fprintf('Final active neurons: %d/%d (%.1f%%)\n', ...
    winner_diversity(epoch), n_pattern_neurons, 100*winner_diversity(epoch)/n_pattern_neurons);
fprintf('Final entropy: %.2f (max possible: %.2f)\n', ...
    winner_entropy(epoch), log(n_pattern_neurons));
fprintf('Total weight change: %.4f\n', norm(W_MT - W_MT_initial, 'fro'));

% Test learned neurons
fprintf('\nTesting pattern selectivity...\n');
test_vels = [-3, 0, 3];
test_dirs = [0, pi/4, pi/2];

response_matrix = zeros(length(test_vels), length(test_dirs), n_pattern_neurons);

for v_idx = 1:length(test_vels)
    for d_idx = 1:length(test_dirs)
        test_stim = create_moving_RDS(test_vels(v_idx), test_dirs(d_idx), 0.7, param);
        C1_test = popFlowV1MT(test_stim.image, param);
        act_test = extract_component_activities_fixed(C1_test, param);
        resp_test = W_MT * act_test;
        
        response_matrix(v_idx, d_idx, :) = resp_test;
        
        [~, best_neurons] = maxk(resp_test, 3);
        fprintf('  V=%.1f, θ=%.2f → Neurons: [%s]\n', ...
            test_vels(v_idx), test_dirs(d_idx), num2str(best_neurons'));
    end
end

%% VISUALIZATION
figure('Position', [100, 100, 1400, 900]);

% Learning curves
subplot(2,4,1);
plot(1:epoch, weight_changes(1:epoch), 'LineWidth', 2);
xlabel('Epoch'); ylabel('Avg Weight Change');
title('Weight Dynamics'); grid on;

subplot(2,4,2);
plot(1:epoch, winner_diversity(1:epoch), 'LineWidth', 2);
xlabel('Epoch'); ylabel('Active Neurons');
title('Population Diversity'); grid on;
ylim([0, n_pattern_neurons]);

subplot(2,4,3);
plot(1:epoch, winner_entropy(1:epoch), 'LineWidth', 2);
hold on;
plot([1, epoch], [log(n_pattern_neurons), log(n_pattern_neurons)], 'r--');
xlabel('Epoch'); ylabel('Entropy');
title('Winner Distribution Entropy'); grid on;
legend('Actual', 'Maximum');

subplot(2,4,4);
bar(winners_count);
xlabel('Pattern Neuron'); ylabel('Times Selected');
title('Final Winner Statistics'); grid on;

% Weight matrices
subplot(2,4,5);
imagesc(W_MT_initial);
colorbar; title('Initial Weights');
xlabel('Component'); ylabel('Pattern Neuron');

subplot(2,4,6);
imagesc(W_MT);
colorbar; title('Final Weights');
xlabel('Component'); ylabel('Pattern Neuron');

subplot(2,4,7);
imagesc(W_MT - W_MT_initial);
colorbar; title('Weight Change');
xlabel('Component'); ylabel('Pattern Neuron');

% Response matrix for one example
subplot(2,4,8);
example_resp = squeeze(response_matrix(:, :, 1));
imagesc(test_dirs, test_vels, example_resp);
colorbar;
xlabel('Direction (rad)'); ylabel('Velocity (px/frame)');
title('Example Neuron Tuning');

% Individual neuron tuning
figure('Position', [100, 100, 1200, 800]);
n_plot = min(12, n_pattern_neurons);
for i = 1:n_plot
    subplot(3, 4, i);
    weights_2d = reshape(W_MT(i, :), [n_orient, n_vel]);
    imagesc(param.pref_vel, 1:n_orient, weights_2d);
    colorbar;
    title(sprintf('Neuron %d', i));
    xlabel('Velocity'); ylabel('Orientation');
end

%% SAVE
save_filename = sprintf('hebbian_topk_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
save(save_filename, 'W_MT', 'W_MT_initial', 'learning_params', ...
    'weight_changes', 'winner_diversity', 'winner_entropy', 'param', ...
    'response_matrix');
fprintf('\nSaved: %s\n', save_filename);

%% HELPER
function activities = extract_component_activities_fixed(C1, param)
    [sy, sx, n_orient, n_vel] = size(C1{1});
    combined = (C1{2});
    activities = squeeze(mean(combined, [1, 2]));
    activities = activities(:);
    activities = max(0, activities);  % Only remove negative values
    activities = activities + 0.01 * max(activities);  % Add small offset
    activities = activities / (norm(activities) + eps);
end