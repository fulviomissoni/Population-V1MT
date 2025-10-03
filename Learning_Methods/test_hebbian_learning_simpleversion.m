%% SIMPLE COMPETITIVE LEARNING - MOST ROBUST
% This is the simplest approach that guarantees all neurons learn

clear; close all; clc;
addpath FUNCTIONS

%% SETUP (same as before - abbreviated)
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
eta_init = 0.3;  % Initial learning rate
eta_final = 0.01;
max_epochs = 300;

% Random initialization
W = randn(n_neurons, n_components) * 0.1;
for i = 1:n_neurons
    W(i,:) = W(i,:) / norm(W(i,:));
end

%% GENERATE TRAINING DATA (minimal set)
fprintf('Generating training data...\n');
vel_vals = linspace(min(param.pref_vel)*0.6, max(param.pref_vel)*0.6, 5);
dir_vals = linspace(0, 2*pi, 10);
[VV, DD] = meshgrid(vel_vals, dir_vals);
VV = VV(:); DD = DD(:);
n_train = length(VV);

activities_all = zeros(n_components, n_train);
for i = 1:n_train
    stim = create_moving_RDS(VV(i), DD(i), 0.6, param);
    C1 = popFlowV1MT(stim.image, param);
    activities_all(:, i) = extract_activities_simple(C1);
    if mod(i, 10) == 0
        fprintf('  %d/%d\n', i, n_train);
    end
end

%% SIMPLE COMPETITIVE LEARNING
fprintf('\nTraining with competitive learning...\n');

winner_counts = zeros(n_neurons, 1);

for epoch = 1:max_epochs
    % Decay learning rate
    eta = eta_init * (eta_final/eta_init)^(epoch/max_epochs);
    
    % Shuffle data
    perm = randperm(n_train);
    
    for i = 1:n_train
        x = activities_all(:, perm(i));
        
        % Find winner (closest neuron)
        distances = sum((W - x').^2, 2);
        [~, winner] = min(distances);
        
        winner_counts(winner) = winner_counts(winner) + 1;
        
        % Update ONLY winner: move toward input
        W(winner, :) = W(winner, :) + eta * (x' - W(winner, :));
        
        % Normalize
        W(winner, :) = W(winner, :) / norm(W(winner, :));
    end
    
    if mod(epoch, 50) == 0
        active = sum(winner_counts > 0);
        fprintf('Epoch %d: η=%.4f, Active=%d/%d\n', epoch, eta, active, n_neurons);
    end
end

%% TEST AND VISUALIZE
fprintf('\nTesting learned representations...\n');

test_vels = [-3, 0, 3];
test_dirs = [0, pi/3, 2*pi/3];

for v_idx = 1:length(test_vels)
    for d_idx = 1:length(test_dirs)
        test_stim = create_moving_RDS(test_vels(v_idx), test_dirs(d_idx), 0.7, param);
        C1_test = popFlowV1MT(test_stim.image, param);
        x_test = extract_activities_simple(C1_test);
        
        % Find best matching neurons
        distances = sum((W - x_test').^2, 2);
        [~, sorted] = sort(distances);
        best_3 = sorted(1:3);
        
        fprintf('  V=%.1f, θ=%.2f → Neurons: %d, %d, %d\n', ...
            test_vels(v_idx), test_dirs(d_idx), best_3(1), best_3(2), best_3(3));
    end
end

%% PLOT
figure('Position', [100, 100, 1400, 600]);

subplot(2,4,1);
bar(winner_counts);
xlabel('Neuron'); ylabel('Times Selected');
title('Winner Distribution');

active_pct = 100 * sum(winner_counts > 0) / n_neurons;
fprintf('\nFinal: %d/%d neurons active (%.1f%%)\n', ...
    sum(winner_counts > 0), n_neurons, active_pct);

% Plot neuron tuning
for i = 1:min(7, n_neurons)
    subplot(2, 4, i+1);
    w_2d = reshape(W(i, :), [n_orient, length(param.pref_vel)]);
    imagesc(param.pref_vel, 1:n_orient, w_2d);
    colorbar;
    title(sprintf('Neuron %d (n=%d)', i, winner_counts(i)));
    xlabel('Velocity'); ylabel('Orientation');
end

%% HELPER
function act = extract_activities_simple(C1)
    combined = (C1{1} + C1{2}) / 2;
    act = squeeze(mean(combined, [1, 2]));
    act = act(:);
    act = act / (norm(act) + eps);
end