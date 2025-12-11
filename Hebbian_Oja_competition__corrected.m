%% FIXED HEBBIAN + OJA WITH COMPETITION
% Key fixes:
% 1. Remove batch z-scoring (inputs already normalized)
% 2. Don't center Yc (keep positive firing rates for Hebbian)
% 3. Softer competition parameters
% 4. Better learning rate schedule
% 5. Oja decay computed correctly

clear; close all; clc;
addpath FUNCTIONS
addpath Learning_Methods\

%% -------------------- SETUP --------------------
filter_file = 'FILTERS/New/Gt43B0.0210f0.063.mat';
k0 = 0.063; samples = 127; n_orient = 8; filter_sample = 43;
ft_choice = 'gabor';
w0_i = linspace(-0.45,0.45,11);
v = w0_i./k0;

param.spat_freq = k0; param.n_orient = n_orient; param.pref_vel = v;
param.temp_filt = ft_choice; param.spatial_filt = filter_file;
param.samples = samples; param.filter_sample = filter_sample;
param.norm_param = [1e-17, 1]; param.sigma_pool = 3; param.num_or_ch_pooled = 8;

sigma_orient = 10; x = linspace(1,n_orient,100); x_8 = round(linspace(1,100,n_orient));
w_o = exp(-(x(x_8)-(1:n_orient)').^2./(2.*sigma_orient.^2));
param.orient_weighting = w_o./sum(w_o(:));

%% -------------------- LEARNING PARAMS --------------------
n_neurons = 25;
n_components = length(param.pref_vel) * param.n_orient;

% FIXED: Better learning rate schedule
eta0 = 3e-3;            % Higher initial LR
eta_min = 1e-6;         % Minimum LR
eta = eta0;
decay_rate = 0.7;     % Smooth exponential decay
anneal = 0.7;           % ogni 'anneal_step' riduco il LR
anneal_step = 1000;

B = 100;                % Batch size (smaller = more stochastic)
T = 80000;              % Total iterations

% FIXED: Softer competition
comp_mode = 'divisive';
alpha = 0.6;            % Reduced from 0.7
sigma_div = 1e-2;       % Reduced from 0.08

% Lateral inhibition (alternative)
gamma = 0.3;
C = ones(n_neurons)/n_neurons; 
C(1:n_neurons+1:end) = 0;

%% -------------------- INITIALIZATION --------------------
fprintf('Random initialization of W (N x M)...\n');
W = randn(n_components, n_neurons) * 0.1;  % Slightly larger init
W = normalize_cols(W);
W_initial = W;

%% -------------------- GENERATE TRAINING DATA --------------------
fprintf('Generating training data...\n');
vel_vals = linspace(param.pref_vel(1),param.pref_vel(end),45);  
dir_vals = linspace(0,pi,30)';
dir_vals(end) = [];

[VV, DD] = meshgrid(vel_vals, dir_vals);
VV = VV(:); DD = DD(:);
n_train = length(VV);

fprintf('Total stimuli: %d\n', n_train);

activities_all = zeros(n_components, n_train);
for i = 1:n_train
    stim = create_moving_RDS(VV(i), DD(i), 0.7, param);
    C1 = popFlowV1MT(stim.image, param);
    activities_all(:, i) = extract_component_activities(C1, param);
    if mod(i,100)==0, fprintf('  %d/%d\n', i, n_train); end
end

% Check stimulus diversity
corr_mat = corrcoef(activities_all');
mean_corr = mean(corr_mat(triu(true(size(corr_mat)),1)));
fprintf('Mean stimulus correlation: %.3f\n', mean_corr);

%% -------------------- HEBBIAN + OJA TRAINING --------------------
fprintf('\nStarting Hebbian+Oja training...\n');
winner_counts = zeros(n_neurons,1);
weight_change_trace = zeros(T,1);
selectivity_trace = zeros(T,1);
lr_trace = zeros(T,1);

for t = 1:T
% sample and batch-normalize (per-channel)
idx = randsample(n_train, B);
Xb = activities_all(:, idx);    % N x B
mu = mean(Xb,2); sd = std(Xb,0,2) + 1e-8;
Xb = (Xb - mu) ./ sd;           % centered & scaled

% linear projection (we still compute it)
Y_lin = W' * Xb;                % M x B

% rectified firing
Y = max(Y_lin, 0);              % M x B

% competition -> post-competition firing Yc (nonnegative)
switch comp_mode
    case 'divisive'
        denom = sigma_div + alpha * sum(Y, 1);   % 1 x B
        Yc = Y ./ denom;                         % M x B
    case 'lateral'
        Yc = Y - gamma * (C * Y);
        Yc = max(Yc, 0);
    otherwise
        error('Unknown comp_mode');
end

% center Yc across the batch to remove DC bias (gives signed signal)
Yc_mean = mean(Yc, 2);          % M x 1
Yc_center = Yc - Yc_mean;       % M x B

% winner counts still from Yc if you want (non-centered)
[~, winners] = max(Yc, [], 1);
for ii = 1:B, winner_counts(winners(ii)) = winner_counts(winners(ii)) + 1; end

% Oja update using centered Yc (both Hebbian and stabilizer use Yc_center)
m2 = mean(Yc_center.^2, 2);     % M x 1
dW = (Xb * Yc_center') / B - W .* (ones(size(W,1),1) * m2');  % N x M
W_old = W;
W = W + eta * dW;

% normalize columns
W = normalize_cols(W);

    % stats
   weight_change_trace(t) = norm(W - W_old, 'fro');

    % anneal learning rate
    if mod(t, anneal_step) == 0
        eta = eta * anneal;
        fprintf(' Iter %d: annealed eta -> %.4e\n', t, eta);
    end
    
    % Progress reporting
    if mod(t, round(T/20)) == 0
        active = sum(winner_counts > 0);
        fprintf(' Iter %d/%d: ΔW=%.4e, LR=%.4e, active=%d/%d, sel=%.3f\n', ...
            t, T, weight_change_trace(t), eta, active, n_neurons, selectivity_trace(t));
    end
end

fprintf('\nTraining finished.\n');

%% -------------------- ANALYSIS --------------------
W_final = W;
active_neurons = sum(winner_counts > 0);
fprintf('\nActive neurons: %d/%d (%.1f%%)\n', active_neurons, n_neurons, 100*active_neurons/n_neurons);
fprintf('Total weight change: %.4f\n', norm(W_final - W_initial, 'fro'));

% Test on training data
final_responses = max(W_final' * activities_all, 0);

% Visualize MT tuning curves
n_dir = length(dir_vals);
n_vel = length(vel_vals);
figure('Name', 'MT Tuning Curves');
popresponse_tiled_countour(permute(reshape(final_responses, 5, 5, n_dir, n_vel), [3,4,1,2]), vel_vals);
sgtitle('Learned MT Neuron Tuning Curves'); colormap('hot')

% Visualize learned weights
figure('Name', 'Learned Weights');
popresponse_tiled_countour(permute(reshape(W_final', 5, 5, n_orient, 11), [3,4,1,2]), param.pref_vel);
sgtitle('Learned Weight Structure');colormap('hot')

%% -------------------- DETAILED DIAGNOSTICS --------------------
figure('Position',[100,100,1400,900]);

% Learning curves
subplot(2,4,1);
semilogy(1:T, weight_change_trace, 'LineWidth', 1.5);
xlabel('Iteration'); ylabel('||ΔW||_F'); 
title('Weight Change (log scale)'); grid on;

subplot(2,4,2);
plot(1:T, lr_trace, 'LineWidth', 1.5);
xlabel('Iteration'); ylabel('Learning Rate');
title('LR Schedule'); grid on;

subplot(2,4,3);
plot(find(selectivity_trace>0), selectivity_trace(selectivity_trace>0), 'LineWidth', 1.5);
xlabel('Iteration'); ylabel('Mean Selectivity');
title('Selectivity Evolution'); grid on;

% Winner distribution
subplot(2,4,4);
bar(winner_counts);
xlabel('Neuron'); ylabel('Times Selected'); 
title(sprintf('Winner Histogram (Active: %d)', active_neurons)); 
grid on;

% Weight matrices
subplot(2,4,5);
imagesc(W_initial'); colorbar; 
title('Initial Weights');
xlabel('V1 Component'); ylabel('MT Neuron');

subplot(2,4,6);
imagesc(W_final'); colorbar; 
title('Final Weights');
xlabel('V1 Component'); ylabel('MT Neuron');

% Weight norms
subplot(2,4,7);
histogram(vecnorm(W_final), 30);
xlabel('Weight Norm'); ylabel('Count');
title('Weight Norm Distribution');
grid on;

% Response distribution
subplot(2,4,8);
histogram(final_responses(:), 50);
xlabel('Response'); ylabel('Count');
title('MT Response Distribution');
grid on;

%% Individual neuron weights
n_vel_param = length(param.pref_vel);
figure('Position',[100,100,1200,1000]);
for i = 1:min(25, n_neurons)
    subplot(5,5,i);
    w_2d = reshape(W_final(:, i), [n_orient, n_vel_param]);
    imagesc(param.pref_vel, 1:n_orient, w_2d);
    colorbar;
    title(sprintf('N%d (n=%d)', i, winner_counts(i)), 'FontSize', 8);
    if i > 20, xlabel('Vel'); end
    if mod(i-1,5)==0, ylabel('Or'); end
end
sgtitle('Individual MT Neuron Weights (Orient × Velocity)');

%% Compute final selectivity
fprintf('\nComputing selectivity indices...\n');
selectivity = zeros(n_neurons, 1);
pref_vel = zeros(n_neurons, 1);
pref_orient = zeros(n_neurons, 1);

for neuron = 1:n_neurons
    w_2d = reshape(W_final(:, neuron), [n_orient, n_vel_param]);
    [max_val, max_idx] = max(w_2d(:));
    mean_val = mean(w_2d(:));
    selectivity(neuron) = (max_val - mean_val) / (max_val + eps);
    
    [or_idx, vel_idx] = ind2sub([n_orient, n_vel_param], max_idx);
    pref_orient(neuron) = or_idx;
    pref_vel(neuron) = param.pref_vel(vel_idx);
end

fprintf('Selectivity: min=%.3f, mean=%.3f, max=%.3f\n', ...
    min(selectivity), mean(selectivity), max(selectivity));
fprintf('Highly selective (>0.4): %d/%d\n', ...
    sum(selectivity > 0.4), n_neurons);

% Velocity preference distribution
figure;
histogram(pref_vel(winner_counts>0), 20);
xlabel('Preferred Velocity'); ylabel('# Neurons');
title('Distribution of Velocity Preferences (Active Neurons)');
grid on;

%% SAVE
save_filename = sprintf('hebbian_oja_FIXED_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
save(save_filename, 'W_final', 'W_initial', 'param', 'winner_counts', ...
    'selectivity', 'activities_all', 'pref_vel', 'pref_orient');
fprintf('\nSaved: %s\n', save_filename);

%% -------------------- HELPER FUNCTIONS --------------------
function Wn = normalize_cols(W)
    norms = vecnorm(W);
    Wn = W ./ max(norms, 1e-12);
end

function sel = compute_mean_selectivity(W, n_orient, n_vel)
    [~, n_neurons] = size(W);
    sel_vals = zeros(n_neurons, 1);
    for i = 1:n_neurons
        w_2d = reshape(W(:,i), [n_orient, n_vel]);
        max_val = max(w_2d(:));
        mean_val = mean(w_2d(:));
        sel_vals(i) = (max_val - mean_val) / (max_val + 1e-8);
    end
    sel = mean(sel_vals);
end

function activities = extract_component_activities(C1, param)
    combined = C1{2};
    activities = squeeze(mean(combined, [1, 2]));
    activities = activities(:);
    activities = max(0, activities);
    activities = activities + 0.001 * max(activities);
    if norm(activities) > 0
        activities = activities / norm(activities);
    end
end