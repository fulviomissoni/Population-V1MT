%% SYSTEMATIC DEBUGGING WORKFLOW FOR HEBBIAN LEARNING
% This script provides diagnostic checks to identify problems

clear
close all
clc
addpath FUNCTIONS

%% LOAD YOUR PARAMETERS (same as before)
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

%% DIAGNOSTIC 1: CHECK C1 RESPONSES ARE NON-ZERO AND VARIED
fprintf('=== DIAGNOSTIC 1: Check C1 Responses ===\n');

% Create a simple test stimulus
test_vel = 2.0;
test_dir = pi/4;
test_coherence = 0.7;

stim = create_moving_RDS(test_vel, test_dir, test_coherence, param);
C1 = popFlowV1MT(stim.image, param);

fprintf('C1{1} size: [%s]\n', num2str(size(C1{1})));
fprintf('C1{2} size: [%s]\n', num2str(size(C1{2})));
fprintf('C1{1} range: [%.4f, %.4f]\n', min(C1{1}(:)), max(C1{1}(:)));
fprintf('C1{2} range: [%.4f, %.4f]\n', min(C1{2}(:)), max(C1{2}(:)));
fprintf('C1{1} mean: %.4f, std: %.4f\n', mean(C1{1}(:)), std(C1{1}(:)));
fprintf('C1{2} mean: %.4f, std: %.4f\n', mean(C1{2}(:)), std(C1{2}(:)));

% Check if responses are all zero or constant
if max(C1{1}(:)) < 1e-10 && max(C1{2}(:)) < 1e-10
    error('PROBLEM: C1 responses are essentially zero! Check stimulus generation.');
end

if std(C1{1}(:)) < 1e-10 || std(C1{2}(:)) < 1e-10
    warning('PROBLEM: C1 responses have no variance! All neurons responding identically.');
end

% Visualize C1 responses
figure('Position', [100, 100, 1200, 500]);
subplot(2,3,1);
imagesc(squeeze(mean(C1{1}, [1,2])));
colorbar; title('C1{1} Tuning (orient x vel)');
xlabel('Velocity Channel'); ylabel('Orientation Channel');

subplot(2,3,2);
imagesc(squeeze(mean(C1{2}, [1,2])));
colorbar; title('C1{2} Tuning (orient x vel)');
xlabel('Velocity Channel'); ylabel('Orientation Channel');

%% DIAGNOSTIC 2: CHECK COMPONENT ACTIVITIES EXTRACTION
fprintf('\n=== DIAGNOSTIC 2: Check Component Activities ===\n');

activities = extract_component_activities(C1, param);
fprintf('Component activities size: %d\n', length(activities));
fprintf('Activities range: [%.4f, %.4f]\n', min(activities), max(activities));
fprintf('Activities mean: %.4f, std: %.4f\n', mean(activities), std(activities));
fprintf('Activities norm: %.4f\n', norm(activities));

% Check if activities are meaningful
if std(activities) < 1e-6
    error('PROBLEM: Component activities have no variance! All values are identical.');
end

if norm(activities) < 1e-10
    error('PROBLEM: Component activities norm is zero!');
end

% Visualize component activities
subplot(2,3,3);
plot(activities, 'o-', 'LineWidth', 2);
title('Component Activities');
xlabel('Component Index'); ylabel('Activity');
grid on;

% Reshape and visualize as tuning
n_vel = length(param.pref_vel);
activities_2d = reshape(activities, [n_orient, n_vel]);
subplot(2,3,4);
imagesc(activities_2d);
colorbar; title('Activities (orient x vel)');
xlabel('Velocity'); ylabel('Orientation');

%% DIAGNOSTIC 3: CHECK WEIGHT INITIALIZATION
fprintf('\n=== DIAGNOSTIC 3: Check Weight Initialization ===\n');

n_pattern_neurons = 88;
n_component_neurons = length(activities);

W_MT_init = initialize_pattern_weights(n_pattern_neurons, n_component_neurons, param);
fprintf('W_MT size: [%d x %d]\n', size(W_MT_init, 1), size(W_MT_init, 2));
fprintf('W_MT range: [%.4f, %.4f]\n', min(W_MT_init(:)), max(W_MT_init(:)));
fprintf('W_MT mean: %.4f, std: %.4f\n', mean(W_MT_init(:)), std(W_MT_init(:)));

% Check weight diversity
row_norms = sqrt(sum(W_MT_init.^2, 2));
fprintf('Weight row norms - mean: %.4f, std: %.4f\n', mean(row_norms), std(row_norms));

if std(W_MT_init(:)) < 1e-6
    error('PROBLEM: All weights are identical after initialization!');
end

% Visualize initial weights
subplot(2,3,5);
imagesc(W_MT_init);
colorbar; title('Initial Weights');
xlabel('Component Neuron'); ylabel('Pattern Neuron');

%% DIAGNOSTIC 4: CHECK PATTERN RESPONSES
fprintf('\n=== DIAGNOSTIC 4: Check Pattern Neuron Responses ===\n');

pattern_responses = W_MT_init * activities;
fprintf('Pattern responses size: %d\n', length(pattern_responses));
fprintf('Pattern responses range: [%.4f, %.4f]\n', min(pattern_responses), max(pattern_responses));
fprintf('Pattern responses mean: %.4f, std: %.4f\n', mean(pattern_responses), std(pattern_responses));

if std(pattern_responses) < 1e-6
    error('PROBLEM: Pattern responses have no variance! Check weight-activity interaction.');
end

% After softmax
pattern_responses_soft = softmax(pattern_responses);
fprintf('After softmax range: [%.4f, %.4f]\n', min(pattern_responses_soft), max(pattern_responses_soft));
fprintf('Softmax sum: %.4f (should be 1.0)\n', sum(pattern_responses_soft));
fprintf('Max softmax value: %.4f (winner strength)\n', max(pattern_responses_soft));

% Visualize pattern responses
subplot(2,3,6);
bar(pattern_responses_soft);
title('Pattern Responses (softmax)');
xlabel('Pattern Neuron'); ylabel('Response');
[~, winner] = max(pattern_responses_soft);
hold on;
plot(winner, pattern_responses_soft(winner), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
text(winner, pattern_responses_soft(winner)*1.1, 'Winner', 'HorizontalAlignment', 'center');

%% DIAGNOSTIC 5: SIMULATE ONE WEIGHT UPDATE
fprintf('\n=== DIAGNOSTIC 5: Simulate One Weight Update ===\n');

W_MT_test = W_MT_init;
[~, winner] = max(pattern_responses_soft);

eta = 1e-3;
weight_decay = 1e-4;

fprintf('Winner neuron: %d\n', winner);
fprintf('Winner weight before update:\n');
fprintf('  First 10 weights: [%s]\n', num2str(W_MT_test(winner, 1:10), '%.4f '));
fprintf('  Weight norm: %.4f\n', norm(W_MT_test(winner, :)));

% Apply Hebbian update
W_before = W_MT_test(winner, :);

for j = 1:n_component_neurons
    pre_activity = activities(j);
    post_activity = pattern_responses_soft(winner);
    
    hebbian_term = pre_activity * post_activity;
    dW = eta * hebbian_term - weight_decay * W_MT_test(winner, j);
    
    W_MT_test(winner, j) = W_MT_test(winner, j) + dW;
end

% Normalize
W_MT_test(winner, :) = W_MT_test(winner, :) / (norm(W_MT_test(winner, :)) + eps);

W_after = W_MT_test(winner, :);

fprintf('Winner weight after update:\n');
fprintf('  First 10 weights: [%s]\n', num2str(W_after(1:10), '%.4f '));
fprintf('  Weight norm: %.4f\n', norm(W_after));
fprintf('  Weight change magnitude: %.6e\n', norm(W_after - W_before));
fprintf('  Max absolute change: %.6e\n', max(abs(W_after - W_before)));

if max(abs(W_after - W_before)) < 1e-10
    error('PROBLEM: Weights are not changing! Weight updates are too small.');
end

%% DIAGNOSTIC 6: CHECK LEARNING RATE SCALING
fprintf('\n=== DIAGNOSTIC 6: Check Learning Rate Scaling ===\n');

% Compute typical Hebbian term magnitude
typical_hebbian = mean(activities) * max(pattern_responses_soft);
fprintf('Typical Hebbian term: %.6e\n', typical_hebbian);
fprintf('Typical weight update (eta * hebbian): %.6e\n', eta * typical_hebbian);
fprintf('Typical weight decay term: %.6e\n', weight_decay * mean(W_MT_init(winner, :)));

% Suggest appropriate learning rate
suggested_eta = 0.1 / typical_hebbian;
fprintf('\nSuggested learning rate: %.6e\n', suggested_eta);

if eta * typical_hebbian < 1e-8
    warning('PROBLEM: Learning rate is too small relative to Hebbian terms!');
    fprintf('Consider increasing eta to at least %.6e\n', suggested_eta);
end

%% DIAGNOSTIC 7: CHECK MULTIPLE STIMULI DIVERSITY
fprintf('\n=== DIAGNOSTIC 7: Check Training Data Diversity ===\n');

n_test_stim = 10;
all_activities = zeros(n_component_neurons, n_test_stim);

for i = 1:n_test_stim
    vel = param.pref_vel(randi(length(param.pref_vel)));
    dir = rand * 2 * pi;
    coh = 0.3 + 0.4 * rand;
    
    stim_temp = create_moving_RDS(vel, dir, coh, param);
    C1_temp = popFlowV1MT(stim_temp.image, param);
    all_activities(:, i) = extract_component_activities(C1_temp, param);
end

% Check diversity across stimuli
activity_means = mean(all_activities, 1);
activity_stds = std(all_activities, 0, 1);

fprintf('Activity means across stimuli: [%.4f - %.4f]\n', min(activity_means), max(activity_means));
fprintf('Activity stds across stimuli: [%.4f - %.4f]\n', min(activity_stds), max(activity_stds));

% Correlation between different stimuli
corr_matrix = corrcoef(all_activities);
mean_corr = mean(corr_matrix(triu(true(size(corr_matrix)), 1)));
fprintf('Mean correlation between different stimuli: %.4f\n', mean_corr);

if mean_corr > 0.9
    warning('PROBLEM: Different stimuli produce highly correlated responses! Need more diversity.');
end

% Visualize stimulus diversity
figure('Position', [100, 100, 1000, 400]);
subplot(1,2,1);
imagesc(all_activities);
colorbar;
title('Component Activities Across Stimuli');
xlabel('Stimulus Index'); ylabel('Component Index');

subplot(1,2,2);
imagesc(corr_matrix);
colorbar;
title('Correlation Matrix Between Stimuli');
xlabel('Stimulus Index'); ylabel('Stimulus Index');

%% DIAGNOSTIC SUMMARY
fprintf('\n=== DIAGNOSTIC SUMMARY ===\n');
fprintf('1. C1 responses: ');
if max(C1{1}(:)) > 1e-10
    fprintf('✓ Non-zero\n');
else
    fprintf('✗ PROBLEM: Zero responses\n');
end

fprintf('2. Component activities: ');
if std(activities) > 1e-6
    fprintf('✓ Diverse\n');
else
    fprintf('✗ PROBLEM: No variance\n');
end

fprintf('3. Weight initialization: ');
if std(W_MT_init(:)) > 1e-6
    fprintf('✓ Diverse\n');
else
    fprintf('✗ PROBLEM: Uniform weights\n');
end

fprintf('4. Pattern responses: ');
if std(pattern_responses) > 1e-6
    fprintf('✓ Diverse\n');
else
    fprintf('✗ PROBLEM: No variance\n');
end

fprintf('5. Weight updates: ');
if max(abs(W_after - W_before)) > 1e-10
    fprintf('✓ Weights changing\n');
else
    fprintf('✗ PROBLEM: Weights not changing\n');
end

fprintf('6. Learning rate: ');
if eta * typical_hebbian > 1e-8
    fprintf('✓ Appropriate scale\n');
else
    fprintf('✗ PROBLEM: Too small\n');
end

fprintf('7. Stimulus diversity: ');
if mean_corr < 0.9
    fprintf('✓ Diverse stimuli\n');
else
    fprintf('✗ PROBLEM: Stimuli too similar\n');
end