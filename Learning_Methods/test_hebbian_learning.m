%% HEBBIAN LEARNING SIMULATION FOR MT PLAID MODEL
% This script implements unsupervised Hebbian learning to derive
% the weighting scheme for MT pattern neurons from component neurons

clear
close all
clc
addpath FUNCTIONS

%% PARAMETERS FROM YOUR ORIGINAL SCRIPT
% Spatial filters
filter_file = 'FILTERS/New/Gt43B0.0210f0.063.mat';
k0 = 0.063;
samples = 127;
n_orient = 8;
filter_sample = 43;

% Temporal filter
ft_choice = 'gabor';
v = linspace(-1,1,11)*round(1/2*1/k0);

% Normalization
sigma_pool = 3;
alpha1 = 1e-17;
alpha2 = 1;

% Organize parameters
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
param.diff_c = linspace(.25,.75,11);

% Orientation weighting
sigma_orient = 10;
x = linspace(1,n_orient,100);
x_8 = round(linspace(1,100,n_orient));
w_o = exp(-(x(x_8)-(1:n_orient)').^2./(2.*sigma_orient.^2));
param.orient_weighting = w_o./sum(w_o(:));
param.sigma_orient = sigma_orient;

%% HEBBIAN LEARNING PARAMETERS
learning_params = struct();
learning_params.eta = 1e-1;             % small learning rate
learning_params.max_epochs = 3000;      % Maximum training epochs
learning_params.conv_threshold = 1e-6;  % Convergence threshold
learning_params.decay_rate = 0.9995;    % Learning rate decay
learning_params.weight_decay = 1e-4;    % L2
learning_params.sigma_w = 0.15;          % target weighting (if used)
% learning_params.batch_size = 32;        % use minibatches if possible



% Initialize MT pattern neurons
n_pattern_neurons = 88;  % Number of pattern neurons to learn
n_component_neurons = length(param.pref_vel) * param.n_orient;

% Better weight initialization (Gaussian around preferred patterns)
W_MT = initialize_pattern_weights(n_pattern_neurons, n_component_neurons, param);

%% TRAINING STIMULUS GENERATION (RDS instead of plaids)
fprintf('Generating RDS training stimuli...\n');

% Create diverse training set with RDS stimuli
[v, orient] = meshgrid(sort(resample(param.pref_vel,40,11)),linspace(0,2*pi,n_orient));
v = v(:); orient = orient(:);
n_train_stimuli = numel(v);
training_data = cell(n_train_stimuli, 1);
target_responses = cell(n_train_stimuli, 1);

for stim_idx = 1:n_train_stimuli
    % Generate random RDS parameters
    true_velocity = v(stim_idx);
    true_direction = orient(stim_idx);  % Full circle
    coherence = 0.001*rand;    % Motion coherence
    
    % Create RDS stimulus
    stim = create_moving_RDS(true_velocity, true_direction, coherence, param);
    
    % Get component responses
    C1 = popFlowV1MT(stim.image, param);
    
    % Store training data
    training_data{stim_idx} = C1;
    target_responses{stim_idx} = [true_velocity, true_direction, coherence];
    
    if mod(stim_idx, 20) == 0
        fprintf('Stimulus %d/%d generated\n', stim_idx, n_train_stimuli);
    end
end

%% HEBBIAN LEARNING ALGORITHM
fprintf('Starting Hebbian learning...\n');

% Learning variables
eta = learning_params.eta;
weight_history = zeros(learning_params.max_epochs, 1);
convergence_history = zeros(learning_params.max_epochs, 1);

for epoch = 1:learning_params.max_epochs
    epoch_error = 0;
    
    % Shuffle training data
    perm = randperm(n_train_stimuli);
    
    for i = 1:n_train_stimuli
        stim_idx = perm(i);
        C1 = training_data{stim_idx};
        true_vel = target_responses{stim_idx}(1);
        true_dir = target_responses{stim_idx}(2);
        
        % Extract component activities (improved function)
        component_activities = extract_component_activities(C1, param);
        
        % Compute pattern neuron responses
        pattern_responses = W_MT * component_activities;
        
        % Apply softmax for competition
        pattern_responses = softmax(pattern_responses);
        
        % Winner-take-all: find most active pattern neuron
        [~, winner_idx] = max(pattern_responses);
        
        % Hebbian update rule (purely unsupervised)
        for j = 1:n_component_neurons
            % Pre-synaptic activity
            pre_activity = component_activities(j);
            
            % Post-synaptic activity (only for winner)
            post_activity = pattern_responses(winner_idx);
            
            % Standard Hebbian term
            hebbian_term = pre_activity * post_activity;
            
            % Weight update with decay
            dW = eta * hebbian_term - learning_params.weight_decay * W_MT(winner_idx, j);
            
            % Update weight
            W_MT(winner_idx, j) = W_MT(winner_idx, j) + dW;
        end
        
        % Normalize weights to prevent runaway growth
        W_MT(winner_idx, :) = W_MT(winner_idx, :) / (norm(W_MT(winner_idx, :)) + eps);
        
        % Compute error for monitoring (optional)
        epoch_error = epoch_error + var(pattern_responses);
    end
    
    % Update learning rate
    eta = eta * learning_params.decay_rate;
    
    % Store convergence metrics
    weight_history(epoch) = mean(abs(W_MT(:)));
    convergence_history(epoch) = epoch_error / n_train_stimuli;
    
    % Print progress
    if mod(epoch, 50) == 0
        fprintf('Epoch %d: Error = %.6f, Mean Weight = %.6f\n', ...
            epoch, convergence_history(epoch), weight_history(epoch));
    end
    
    % Check convergence
    if epoch > 10 && abs(convergence_history(epoch) - convergence_history(epoch-1)) < learning_params.conv_threshold
        fprintf('Converged at epoch %d\n', epoch);
        break;
    end
end

% %% VISUALIZE POPULATION RESPONSES
% fprintf('Visualizing population responses...\n');
% 
% % Test with sample velocities and directions
% test_velocities = linspace(min(param.pref_vel), max(param.pref_vel), 15);
% test_directions = linspace(0, 2*pi, 16);
% test_coherence = 0.7;
% 
% % Initialize response matrix
% population_responses = zeros(length(test_velocities), length(test_directions), n_pattern_neurons);
% 
% for v_idx = 1:length(test_velocities)
%     for d_idx = 1:length(test_directions)
%         % Create test stimulus
%         test_vel = test_velocities(v_idx);
%         test_dir = test_directions(d_idx);
% 
%         stim = create_moving_RDS(test_vel, test_dir, test_coherence, param);
%         C1_test = popFlowV1MT(stim.image, param);
% 
%         % Extract component activities
%         component_activities_test = extract_component_activities(C1_test, param);
% 
%         % Compute pattern neuron responses
%         pattern_responses_test = W_MT * component_activities_test;
% 
%         % Store responses
%         population_responses(v_idx, d_idx, :) = pattern_responses_test;
%     end
% end
% 
% % Create surf plot for each pattern neuron
% plot_pattern_neuron_responses(population_responses, test_velocities, test_directions, param);
% 
% % Also plot average population response
% plot_average_population_response(population_responses, test_velocities, test_directions);

%% ANALYSIS AND VISUALIZATION
fprintf('Analyzing learned weights...\n');

% Plot learning curves
figure('Position', [100, 100, 1200, 400]);
subplot(1, 3, 1);
plot(1:epoch, convergence_history(1:epoch), 'LineWidth', 2);
xlabel('Epoch'); ylabel('RMS Error');
title('Convergence History');
grid on;

subplot(1, 3, 2);
plot(1:epoch, weight_history(1:epoch), 'LineWidth', 2);
xlabel('Epoch'); ylabel('Mean Weight Magnitude');
title('Weight Evolution');
grid on;

% Analyze weight distributions
subplot(1, 3, 3);
histogram(W_MT(:), 50);
xlabel('Weight Value'); ylabel('Frequency');
title('Learned Weight Distribution');
grid on;

%% SAVE RESULTS
save_filename = sprintf('hebbian_results_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
save(save_filename, 'W_MT', 'learning_params', 'convergence_history', 'weight_history', 'param', 'training_data');
fprintf('Results saved to: %s\n', save_filename);