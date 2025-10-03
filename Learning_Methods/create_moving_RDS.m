function stim = create_moving_RDS(velocity, direction, coherence, param)
    % Create proper moving Random Dot Stereogram stimulus
    % Returns monocular random dot motion stimulus
    
    frames = 43;
    [rows, cols] = deal(param.samples, param.samples);
    
    % Initialize with zeros (black background)
    image_seq = zeros(rows, cols, frames);
    
    % Dot parameters
    dot_density = 0.2;      % Increased density for better stimulation
    dot_size = 2;           % Smaller dots work better
    n_dots = round(rows * cols * dot_density);
    
    % Initialize dot positions (random across space)
    dot_positions = [randi(rows, n_dots, 1), randi(cols, n_dots, 1)];
    
    % Determine coherent vs incoherent dots
    n_coherent = round(n_dots * coherence);
    coherent_mask = false(n_dots, 1);
    coherent_mask(1:n_coherent) = true;
    
    % Calculate displacement per frame
    shift_x = velocity * cos(direction);
    shift_y = velocity * sin(direction);
    
    for f = 1:frames
        frame = zeros(rows, cols);
        
        if f > 1
            % Move coherent dots
            dot_positions(coherent_mask, 1) = dot_positions(coherent_mask, 1) + shift_y;
            dot_positions(coherent_mask, 2) = dot_positions(coherent_mask, 2) + shift_x;
            
            % Move incoherent dots randomly
            random_shifts = randn(sum(~coherent_mask), 2) * abs(velocity);
            dot_positions(~coherent_mask, :) = dot_positions(~coherent_mask, :) + random_shifts;
            
            % Wrap dots around screen boundaries
            dot_positions(:,1) = mod(dot_positions(:,1)-1, rows) + 1;
            dot_positions(:,2) = mod(dot_positions(:,2)-1, cols) + 1;
            
            % Randomly regenerate some incoherent dots (dynamic noise)
            regenerate_idx = find(~coherent_mask);
            regenerate_idx = regenerate_idx(rand(size(regenerate_idx)) < 0.1);
            dot_positions(regenerate_idx, :) = [randi(rows, length(regenerate_idx), 1), ...
                                              randi(cols, length(regenerate_idx), 1)];
        end
        
        % Efficient dot rendering using linear indices
        y_pos = round(dot_positions(:,1));
        x_pos = round(dot_positions(:,2));
        
        % Keep dots within bounds
        valid = y_pos >= 1 & y_pos <= rows & x_pos >= 1 & x_pos <= cols;
        y_pos = y_pos(valid);
        x_pos = x_pos(valid);
        
        % Set dot pixels to 1 (white)
        linear_indices = sub2ind([rows, cols], y_pos, x_pos);
        frame(linear_indices) = 1;
        
        % Optional: Add anti-aliasing for smoother motion
        % frame = imgaussfilt(frame, 0.5);
        
        image_seq(:, :, f) = frame;
    end
    
    stim.image = image_seq;
    stim.velocity = velocity;
    stim.direction = direction;
    stim.coherence = coherence;
    stim.description = 'Monocular RDS with coherent motion';
end