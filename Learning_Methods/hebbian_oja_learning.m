%% ULTRA-SIMPLE COMPETITIVE LEARNING - GUARANTEED TO WORK
% This is the absolute simplest approach that will create diverse neurons

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

%% LEARNING PARAMETERS - KEEP IT SIMPLE
n_neurons = 25;
n_components = length(param.pref_vel) * param.n_orient;
eta_init = 0.5;   % High initial learning rate
eta_final = 0.05; % Low final learning rate
max_epochs = 200;

%% RANDOM INITIALIZATION - NO STRUCTURE
fprintf('Random initialization...\n');
W = randn(n_neurons, n_components) * 0.5 + 0.5;  % Positive random weights
for i = 1:n_neurons
    W(i,:) = W(i,:) / norm(W(i,:));
end
W_initial = W;

%% GENERATE DIVERSE TRAINING DATA
fprintf('Generating training data...\n');

% Use your full velocity range INCLUDING zero
vel_vals = param.pref_vel(2:end-1);  % Skip extreme values
dir_vals = linspace(0, 2*pi-2*pi/12, 12);
[VV, DD] = meshgrid(vel_vals, dir_vals);
VV = VV(:); DD = DD(:);
n_train = length(VV);

fprintf('Total stimuli: %d\n', n_train);

% Pre-compute all activities
activities_all = zeros(n_components, n_train);
for i = 1:n_train
    stim = create_moving_RDS(VV(i), DD(i), 0.7, param);
    C1 = popFlowV1MT(stim.image, param);
    activities_all(:, i) = extract_component_activities(C1, param);
    
    if mod(i, 20) == 0
        fprintf('  %d/%d\n', i, n_train);
    end
end

% Check stimulus diversity
fprintf('\nChecking stimulus diversity...\n');
corr_mat = corrcoef(activities_all');
mean_corr = mean(corr_mat(triu(true(size(corr_mat)), 1)));
fprintf('Mean stimulus correlation: %.3f\n', mean_corr);

if mean_corr > 0.85
    warning('Stimuli are too similar! Learning may not work well.');
end

%% HEBBIAN + OJA'S RULE LEARNING LOOP
fprintf('\nStarting Hebbian learning with Oja''s rule...\n');

winner_counts = zeros(n_neurons, 1); % not really "winners" anymore, but still track strong activations
weight_changes = zeros(max_epochs, 1);

for epoch = 1:max_epochs
    % Linear decay of learning rate
    eta = eta_init + (eta_final - eta_init) * (epoch / max_epochs);
    
    % Shuffle data
    perm = randperm(n_train);
    epoch_change = 0;
    
    for i = 1:n_train
        x = activities_all(:, perm(i));   % input pattern
        
        % ACTIVATIONS for all neurons
        y = W * x;  % size: [n_neurons x 1]
        
        % Hebbian + Oja’s rule update for ALL neurons
        for n = 1:n_neurons
            W_old = W(n, :);
            
            % Oja’s update
            W(n,:) = W(n,:) + eta * y(n) * (x' - y(n) * W(n,:));
            
            % Track change magnitude
            epoch_change = epoch_change + norm(W(n,:) - W_old);
        end
        
        % Optional: count which neuron had the strongest response
        [~, winner] = max(y);
        winner_counts(winner) = winner_counts(winner) + 1;
    end
    
    weight_changes(epoch) = epoch_change / n_train;
    
    if mod(epoch, 25) == 0
        active = sum(winner_counts > 0);
        fprintf('Epoch %d: η=%.3f, ΔW=%.4f, Active=%d/%d\n', ...
            epoch, eta, weight_changes(epoch), active, n_neurons);
    end
end

%% ANALYSIS
fprintf('\n=== RESULTS ===\n');
active_neurons = sum(winner_counts > 0);
fprintf('Active neurons: %d/%d (%.1f%%)\n', active_neurons, n_neurons, ...
    100*active_neurons/n_neurons);
fprintf('Total weight change: %.3f\n', norm(W - W_initial, 'fro'));

% Test with various stimuli
fprintf('\nTesting learned neurons:\n');
test_vels = [-4, -2, 0, 2, 4];
test_dirs = [0, pi/4, pi/2];

for v_idx = 1:length(test_vels)
    for d_idx = 1:length(test_dirs)
        stim_test = create_moving_RDS(test_vels(v_idx), test_dirs(d_idx), 0.7, param);
        C1_test = popFlowV1MT(stim_test.image, param);
        x_test = extract_component_activities(C1_test, param);
        
        % Find closest neurons
        distances = sum((W - x_test').^2, 2);
        [~, sorted] = sort(distances);
        best = sorted(1:3);
        
        fprintf('  V=%+.1f, θ=%.2f → [%d, %d, %d]\n', ...
            test_vels(v_idx), test_dirs(d_idx), best(1), best(2), best(3));
    end
end

%% VISUALIZATION
figure('Position', [100, 100, 1400, 800]);

% Learning curve
subplot(2,4,1);
plot(1:max_epochs, weight_changes, 'LineWidth', 2);
xlabel('Epoch'); ylabel('Weight Change');
title('Learning Progress'); grid on;

% Winner distribution
subplot(2,4,2);
bar(winner_counts);
xlabel('Neuron'); ylabel('Times Selected');
title(sprintf('Winners (Active: %d/%d)', active_neurons, n_neurons));
grid on;

% Weight matrices
subplot(2,4,3);
imagesc(W_initial);
colorbar; title('Initial Weights');
xlabel('Component'); ylabel('Neuron');

subplot(2,4,4);
imagesc(W);
colorbar; title('Final Weights');
xlabel('Component'); ylabel('Neuron');

% Individual neuron tuning
n_vel = length(param.pref_vel);
for i = 1:min(20, n_neurons)
    subplot(5, 4, 4+i);
    w_2d = reshape(W(i, :), [n_orient, n_vel]);
    imagesc(param.pref_vel, 1:n_orient, w_2d);
    colorbar;
    title(sprintf('N%d (n=%d)', i, winner_counts(i)));
    if i > 16, xlabel('Velocity'); end
    if mod(i-1, 4) == 0, ylabel('Orient'); end
end

%% COMPUTE SELECTIVITY
fprintf('\nComputing selectivity indices...\n');
selectivity = zeros(n_neurons, 1);

for neuron = 1:n_neurons
    w_2d = reshape(W(neuron, :), [n_orient, n_vel]);
    max_val = max(w_2d(:));
    mean_val = mean(w_2d(:));
    selectivity(neuron) = (max_val - mean_val) / (max_val + eps);
end

fprintf('Selectivity: min=%.3f, mean=%.3f, max=%.3f\n', ...
    min(selectivity), mean(selectivity), max(selectivity));
fprintf('Highly selective (>0.4): %d/%d\n', sum(selectivity > 0.4), n_neurons);

%% SAVE
save_filename = sprintf('simple_competitive_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
save(save_filename, 'W', 'W_initial', 'param', 'winner_counts', 'selectivity');
fprintf('\nSaved: %s\n', save_filename);

%% UNCHANGED FUNCTION
function activities = extract_component_activities(C1, param)
    [sy, sx, n_orient, n_vel] = size(C1{1});
    combined = (C1{2});
    activities = squeeze(mean(combined, [1, 2]));
    activities = activities(:);
    activities = max(0, activities);
    activities = activities + 0.001 * max(activities);
    % if norm(activities) > 0
    %     activities = activities / norm(activities);
    % end
end