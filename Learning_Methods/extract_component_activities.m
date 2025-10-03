%%
function activities = extract_component_activities(C1, param)
    % Extract and properly format component activities from C1 responses
    
    [sy, sx, n_orient, n_vel] = size(C1{1});
    
    % Average over spatial position and time to get tuning curves
    activities = zeros(n_orient, n_vel);
    
    for i = 1:2  % Both energy channels
        % Average over space and time
        response_channel = squeeze(mean(C1{i}, [1,2]));
        activities = activities + squeeze(response_channel);
    end
    
    % Flatten to vector: order is velocity first, then orientation
    activities = activities(:);
    
    % Normalize to unit length
    activities = activities / (norm(activities) + eps);
end