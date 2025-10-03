function plot_pattern_neuron_responses(population_responses, test_velocities, test_directions, param)
    % Plot surf plots for pattern neuron responses
    
    [n_vel, n_dir, n_neurons] = size(population_responses);
    
    % Plot first 9 neurons in a 3x3 grid
    figure('Position', [100, 100, 1200, 900]);
    sgtitle('Pattern Neuron Tuning Surfaces', 'FontSize', 16);
    
    n_plot = min(9, n_neurons);
    plot_idx = 1;
    
    for i = 1:n_plot
        subplot(3, 3, plot_idx);
        
        % Convert to Cartesian coordinates for better visualization
        [DIR, VEL] = meshgrid(test_directions, test_velocities);
        [X, Y] = pol2cart(DIR, VEL);
        
        responses = squeeze(population_responses(:, :, i));
        
        surf(X, Y, responses, 'EdgeColor', 'none');
        xlabel('Vx (pixels/frame)');
        ylabel('Vy (pixels/frame)');
        zlabel('Response');
        title(sprintf('Neuron %d', i));
        colorbar;
        
        plot_idx = plot_idx + 1;
    end
end