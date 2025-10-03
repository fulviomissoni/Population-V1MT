function target_weights = compute_target_weights(component_activities, true_velocity, true_direction, param)
    % Compute target weights using your original weighting formula
    % W_MT = exp(-W_MT^2/2*sigma^2) where sigma = 0.25
    % W_MT = V_true * [Theta_V1_MT - Theta_true] - V_V1_MT
    
    sigma = 0.25;
    n_components = length(component_activities);
    target_weights = zeros(n_components, 1);
    
    % Get preferred velocities and orientations for each component
    n_vel = length(param.pref_vel);
    n_orient = param.n_orient;
    
    for i = 1:n_components
        % Determine component's preferred velocity and orientation
        vel_idx = ceil(i / n_orient);
        orient_idx = mod(i-1, n_orient) + 1;
        
        V_V1_MT = param.pref_vel(vel_idx);
        Theta_V1_MT = (orient_idx - 1) * pi / n_orient;  % Convert to radians
        
        % Compute weight using your formula
        W_MT_val = true_velocity * (Theta_V1_MT - true_direction) - V_V1_MT;
        target_weights(i) = exp(-W_MT_val^2 / (2 * sigma^2));
    end
    
    % Normalize target weights
    target_weights = target_weights / (sum(target_weights) + eps);
end