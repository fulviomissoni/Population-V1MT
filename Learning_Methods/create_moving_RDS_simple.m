function stim = create_moving_RDS_simple(velocity, direction, coherence, param)
    % Ultra-simple version for debugging
    
    frames = 20;  % Fewer frames for speed
    [rows, cols] = deal(param.samples, param.samples);
    
    image_seq = zeros(rows, cols, frames);
    
    % Create drifting grating as simple alternative
    [x, y] = meshgrid(1:cols, 1:rows);
    
    for f = 1:frames
        % Moving sine wave grating
        phase = 2*pi*(x*cos(direction) + y*sin(direction)) / 20 + f*velocity*0.5;
        frame = 0.5 + 0.5*sin(phase);
        
        % Add noise for robustness
        frame = frame + randn(rows, cols) * 0.1;
        frame = max(0, min(1, frame));  % Clip to [0,1]
        
        image_seq(:, :, f) = frame;
    end
    
    stim.image = image_seq;
    stim.velocity = velocity;
    stim.direction = direction;
    stim.coherence = coherence;
end