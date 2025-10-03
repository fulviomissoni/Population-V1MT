function stim = create_plaid_stimulus(velocity, direction, contrast, param)
    % Create a simple plaid stimulus for training
    stim = struct();
    
    % Generate random dot stimulus moving in specified direction
    frames = 43;
    [rows, cols] = deal(param.samples, param.samples);
    
    % Initialize image sequence
    image_seq = zeros(rows, cols, frames);
    
    % Generate random dots
    dot_density = 0.3;
    for f = 1:frames
        % Create random dot pattern
        dots = rand(rows, cols) < dot_density;
        
        % Apply motion by shifting pattern
        shift_x = round(velocity * cos(direction) * (f-1));
        shift_y = round(velocity * sin(direction) * (f-1));
        
        % Circular shift to simulate motion
        dots_shifted = circshift(dots, [shift_y, shift_x]);
        
        % Apply contrast
        image_seq(:, :, f) = dots_shifted * contrast;
    end
    
    stim.image = image_seq;
    stim.velocity = velocity;
    stim.direction = direction;
    stim.contrast = contrast;
end