function W_MT = initialize_pattern_weights(n_pattern_neurons, n_component_neurons, param)
    % Initialize weights with biologically plausible structure
    
    n_vel = length(param.pref_vel);
    n_orient = param.n_orient;
    
    W_MT = zeros(n_pattern_neurons, n_component_neurons);
    
    for i = 1:n_pattern_neurons
        % Assign each pattern neuron a preferred velocity and direction
        pref_vel_idx = randi(n_vel);
        pref_orient_idx = randi(n_orient);
        
        % Create Gaussian tuning around preferred feature
        for vel_idx = 1:n_vel
            for orient_idx = 1:n_orient
                comp_idx = (vel_idx-1)*n_orient + orient_idx;
                
                % Distance in feature space
                vel_dist = abs(vel_idx - pref_vel_idx) / n_vel;
                orient_dist = min([abs(orient_idx - pref_orient_idx), ...
                                 n_orient - abs(orient_idx - pref_orient_idx)]) / n_orient;
                
                % Gaussian tuning
                W_MT(i, comp_idx) = exp(-(vel_dist^2 + orient_dist^2) / 2*0.001^2);
            end
        end
        
        % Normalize weights
        W_MT(i, :) = W_MT(i, :) / (norm(W_MT(i, :)) + eps);
    end
end