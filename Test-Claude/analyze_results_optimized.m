function analyze_results_optimized(results, param1_name, param1_vals, param2_name, param2_vals)
% Optimized analysis function with memory-efficient processing

fprintf('Analyzing results...\n');

% Extract response data efficiently
n_conditions = length(results.responses);
response_matrix = zeros(length(param2_vals), length(param1_vals));

for i = 1:n_conditions
    if ~isempty(results.responses{i})
        % Extract relevant metric (e.g., peak response)
        response_data = results.responses{i};
        if iscell(response_data)
            metric = max(response_data{1}(:)); % Example: peak response
        else
            metric = max(response_data(:));
        end
        
        % Map to grid position
        [row, col] = ind2sub(size(response_matrix), i);
        response_matrix(row, col) = metric;
    end
end

% Create visualization
figure('Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
imagesc(param1_vals, param2_vals, response_matrix);
colorbar;
xlabel(param1_name);
ylabel(param2_name);
title('Response Heatmap');

subplot(1, 3, 2);
plot(param1_vals, mean(response_matrix, 1), 'LineWidth', 2);
xlabel(param1_name);
ylabel('Mean Response');
title(['Response vs ' param1_name]);
grid on;

subplot(1, 3, 3);
plot(param2_vals, mean(response_matrix, 2), 'LineWidth', 2);
xlabel(param2_name);
ylabel('Mean Response');
title(['Response vs ' param2_name]);
grid on;

% Save analysis results
save('analysis_results.mat', 'response_matrix', 'param1_vals', 'param2_vals');
fprintf('Analysis completed and saved.\n');
end