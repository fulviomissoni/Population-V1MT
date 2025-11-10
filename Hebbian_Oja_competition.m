%% HEBBIAN + OJA WITH COMPETITION (adapted to your V1->MT pipeline)
% Merges your original dataset generation & visualization with a Hebbian+Oja
% training loop using mini-batches and either divisive or lateral competition.
% Author: adapted for user
clear; 
close all; clc;
addpath FUNCTIONS
addpath Learning_Methods\

%% -------------------- SETUP (from your original script) --------------------
filter_file = 'FILTERS/New/Gt43B0.0210f0.063.mat';
k0 = 0.063; samples = 127; n_orient = 8; filter_sample = 43;
ft_choice = 'gabor';
w0_i = linspace(-0.45,0.45,11);
v = w0_i./k0;
% v = linspace(-1,1,11)*round(1/2*1/k0);

param.spat_freq = k0; param.n_orient = n_orient; param.pref_vel = v;
param.temp_filt = ft_choice; param.spatial_filt = filter_file;
param.samples = samples; param.filter_sample = filter_sample;
param.norm_param = [1e-17, 1]; param.sigma_pool = 3; param.num_or_ch_pooled = 8;

sigma_orient = 10; x = linspace(1,n_orient,100); x_8 = round(linspace(1,100,n_orient));
w_o = exp(-(x(x_8)-(1:n_orient)').^2./(2.*sigma_orient.^2));
param.orient_weighting = w_o./sum(w_o(:));

%% learning params (adjust if you want)
n_neurons = 88;                         % M (MT neurons)
n_components = length(param.pref_vel) * param.n_orient;  % input dimensionality (N)
eta0 = 3e-3;            % initial learning rate
eta = eta0;
anneal = 0.9;           % ogni 'anneal_step' riduco il LR
anneal_step = 1000;

B = 150;                % mini-batch size
T = 100000;               % total mini-batch iterations

% competition mode: 'divisive' or 'lateral'
comp_mode = 'divisive';

% divisive competition params
alpha = 0.7;        % ridotto da 0.6
sigma_div = 8e-2;    % aumentato da 1e-2
% alpha = 0.6;
% sigma_div = 1e-2;

% lateral inhibition params (used if comp_mode = 'lateral')
gamma = 0.5;
C = ones(n_neurons)/n_neurons; C(1:n_neurons+1:end) = 0;

%% -------------------- RANDOM INITIALIZATION --------------------
% NOTE: W is N x M (n_components x n_neurons) to match Oja formulation
fprintf('Random initialization of W (N x M)...\n');
W = randn(n_components, n_neurons) * 0.01;
W = normalize_cols(W);
W_initial = W;

%% -------------------- GENERATE TRAINING DATA (reuse your code) --------------------
fprintf('Generating training data...\n');
vel_vals = linspace(param.pref_vel(1),param.pref_vel(end),42);  
dir_vals = linspace(0,pi,31)';%2*pi*rand(N,1);
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
    if mod(i,20)==0, fprintf('  %d/%d\n', i, n_train); end
end

% stimulus similarity check
corr_mat = corrcoef(activities_all');
mean_corr = mean(corr_mat(triu(true(size(corr_mat)),1)));
fprintf('Mean stimulus correlation: %.3f\n', mean_corr);
if mean_corr > 0.85
    warning('Stimuli are very similar: consider enlarging stimulus set.');
end

%% -------------------- HEBBIAN + OJA TRAINING LOOP --------------------
fprintf('\nStarting Hebbian+Oja training (mini-batches)...\n');
winner_counts = zeros(n_neurons,1);
weight_change_trace = zeros(T,1);

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
    weight_change_trace(t) = norm(W - W_old, 'fro') / (n_comp_val(W)+eps);

    % anneal learning rate
    if mod(t, anneal_step) == 0
        eta = eta * anneal;
        fprintf(' Iter %d: annealed eta -> %.4e\n', t, eta);
    end

    if mod(t, round(T/10)) == 0
        fprintf(' Iter %d/%d: ΔW=%.6f, active=%d/%d\n', t, T, weight_change_trace(t), sum(winner_counts>0), n_neurons);
    end
end

fprintf('\nTraining finished.\n');

%% -------------------- ANALYSIS & VISUALIZATION --------------------
% convert W for plotting compatibility: each neuron -> W(:,m)
W_final = W;
active_neurons = sum(winner_counts > 0);
fprintf('Active neurons after Oja training: %d/%d (%.1f%%)\n', active_neurons, n_neurons, 100*active_neurons/n_neurons);
fprintf('Total weight change (Frobenius): %.4f\n', norm(W_final - W_initial, 'fro'));

final_cells = max(W_final' * activities_all, 0);
% figure, popresponse_tiled(permute(reshape(final_cells,5,5,numel(dir_vals),numel(vel_vals)),[3,4,1,2]),vel_vals); title('MT tuning curve');
figure, popresponse_tiled_countour(permute(reshape(final_cells,8,11,numel(dir_vals),numel(vel_vals)),[3,4,1,2]),vel_vals); title('MT tuning curve');
% levels = linspace(min(W_final(:)), max(W_final(:)), 10);

figure, popresponse_tiled_countour(permute(reshape(W_final',8,11,n_orient,numel(param.pref_vel)),[3,4,1,2]),param.pref_vel); title('Learned Weight');

%%

%%
% Learning curve
figure('Position',[100,100,1200,700]);
subplot(2,4,1);
plot(1:T, weight_change_trace, 'LineWidth', 2);
xlabel('Iteration'); ylabel('Avg ΔW'); title('Learning Progress'); grid on;

% Winner distribution
subplot(2,4,2);
bar(winner_counts);
xlabel('Neuron'); ylabel('Times Selected'); title(sprintf('Winners (Active: %d/%d)', active_neurons, n_neurons)); grid on;

% Weight matrices (initial vs final)
subplot(2,4,3);
imagesc(W_initial'); colorbar; title('Initial Weights (neurons x components)');
xlabel('Component'); ylabel('Neuron');

subplot(2,4,4);
imagesc(W_final'); colorbar; title('Final Weights (neurons x components)');
xlabel('Component'); ylabel('Neuron');

% Individual neuron tuning (reshape into orient x vel)
n_vel = length(param.pref_vel);
for i = 1:min(25, n_neurons)
    subplot(5, 5, i);
    w_2d = reshape(W_final(:, i), [n_orient, n_vel]);
    imagesc(param.pref_vel, 1:n_orient, w_2d);
    colorbar;
    title(sprintf('N%d (n=%d)', i, winner_counts(i)));
    if i > 16, xlabel('Velocity'); end
    if mod(i-1, 4) == 0, ylabel('Orient'); end
end
figure, plot(weight_change_trace)
%% compute selectivity index (same metric as before)
fprintf('\nComputing selectivity indices...\n');
selectivity = zeros(n_neurons, 1);
for neuron = 1:n_neurons
    w_2d = reshape(W_final(:, neuron), [n_orient, n_vel]);
    max_val = max(w_2d(:));
    mean_val = mean(w_2d(:));
    selectivity(neuron) = (max_val - mean_val) / (max_val + eps);
end
fprintf('Selectivity: min=%.3f, mean=%.3f, max=%.3f\n', min(selectivity), mean(selectivity), max(selectivity));
fprintf('Highly selective (>0.4): %d/%d\n', sum(selectivity > 0.4), n_neurons);

%% SAVE
save_filename = sprintf('hebbian_oja_competitive_%s.mat', datetime("now"));
save(save_filename, 'W_final', 'W_initial', 'param', 'winner_counts', 'selectivity','activities_all');
fprintf('Saved: %s\n', save_filename);

%% -------------------- HELPER FUNCTIONS --------------------
function Wn = normalize_cols(W)
    % normalize columns to unit norm
    norms = vecnorm(W);
    Wn = W ./ max(norms, 1e-12);
end

function Xz = zscore_batch(X)
    % z-score across batch for each channel (row)
    mu = mean(X, 2);
    sd = std(X, 0, 2) + 1e-8;
    Xz = (X - mu) ./ sd;
end

function n = n_comp_val(W)
    % small helper: number of components used in normalization denominator
    n = numel(W);
end

%% -------------------- UNCHANGED FUNCTION (from your original file) --------------------
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
