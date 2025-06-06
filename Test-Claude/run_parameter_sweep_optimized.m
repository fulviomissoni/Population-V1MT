function results = run_parameter_sweep_optimized(base_stim, base_param, varargin)
% Optimized parameter sweep for visual pathway simulations
% Uses vectorization and memory-efficient processing

% Parse input parameters
p = inputParser;
addParameter(p, 'contr', [], @isnumeric);
addParameter(p, 'norm_alpha', [], @isnumeric);
addParameter(p, 'batch_size', 50, @isnumeric); % Process in batches to manage memory
addParameter(p, 'save_intermediate', true, @islogical);
addParameter(p, 'output_dir', 'results', @ischar);
parse(p, varargin{:});

contr_vals = p.Results.contr;
alpha_vals = p.Results.norm_alpha;
batch_size = p.Results.batch_size;

% Create parameter grid
[contr_grid, alpha_grid] = meshgrid(contr_vals, alpha_vals);
n_total = numel(contr_grid);

% Initialize results structure
results = struct();
results.contr = contr_grid(:);
results.alpha = alpha_grid(:);
results.responses = cell(n_total, 1);
results.timing = zeros(n_total, 1);

% Create output directory if saving intermediate results
if p.Results.save_intermediate && ~exist(p.Results.output_dir, 'dir')
    mkdir(p.Results.output_dir);
end

% Process in batches to manage memory
n_batches = ceil(n_total / batch_size);
fprintf('Processing %d parameter combinations in %d batches...\n', n_total, n_batches);

for batch_idx = 1:n_batches
    batch_start = (batch_idx - 1) * batch_size + 1;
    batch_end = min(batch_idx * batch_size, n_total);
    batch_indices = batch_start:batch_end;
    
    fprintf('Processing batch %d/%d (indices %d-%d)...\n', ...
        batch_idx, n_batches, batch_start, batch_end);
    
    % Process batch using parfeval for asynchronous execution
    if exist('gcp', 'file') && ~isempty(gcp('nocreate'))
        futures = cell(length(batch_indices), 1);
        
        for i = 1:length(batch_indices)
            idx = batch_indices(i);
            
            % Create stimulus with current parameters
            current_stim = base_stim;
            current_stim.contrast_g = [contr_grid(idx), contr_grid(idx)];
            
            % Create parameter set with current alpha
            current_param = base_param;
            current_param.norm_param(2, :) = alpha_grid(idx);
            
            % Submit job asynchronously
            futures{i} = parfeval(@process_single_condition, 1, ...
                current_param, current_stim, idx);
        end
        
        % Collect results
        for i = 1:length(batch_indices)
            idx = batch_indices(i);
            [completed_idx, response, timing] = fetchNext(futures{:});
            
            results.responses{completed_idx} = response;
            results.timing(completed_idx) = timing;
            
            fprintf('Completed condition %d/%d\n', completed_idx, n_total);
        end
    else
        % Sequential processing if no parallel pool
        for i = 1:length(batch_indices)
            idx = batch_indices(i);
            
            current_stim = base_stim;
            current_stim.contrast_g = [contr_grid(idx), contr_grid(idx)];
            
            current_param = base_param;
            current_param.norm_param(2, :) = alpha_grid(idx);
            
            [response, timing] = process_single_condition(current_param, current_stim, idx);
            
            results.responses{idx} = response;
            results.timing(idx) = timing;
            
            fprintf('Completed condition %d/%d\n', idx, n_total);
        end
    end
    
    % Save intermediate results
    if p.Results.save_intermediate
        batch_filename = fullfile(p.Results.output_dir, ...
            sprintf('batch_%03d.mat', batch_idx));
        batch_results = results;
        batch_results.responses = results.responses(batch_start:batch_end);
        batch_results.timing = results.timing(batch_start:batch_end);
        save(batch_filename, 'batch_results');
    end
    
    % Clear memory
    clear futures;
    if exist('gcp', 'file') && ~isempty(gcp('nocreate'))
        parfevalOnAll(@clearvars);
    end
end

% Save final results
if p.Results.save_intermediate
    save(fullfile(p.Results.output_dir, 'final_results.mat'), 'results');
end

fprintf('Parameter sweep completed. Total time: %.2f seconds\n', sum(results.timing));
end



