function results = run_batched_sweep(II, base_param, param1_name, param1_values, param2_name, param2_values)
% Implementazione con batch processing per casi con memoria limitata

n_param1 = length(param1_values);
n_param2 = length(param2_values);
batch_size = 10; % Processa 10 combinazioni alla volta

results = cell(n_param1, n_param2, 2);
total_combinations = n_param1 * n_param2;

fprintf('Elaborazione a batch (dimensione batch: %d)...\n', batch_size);

h = waitbar(0, 'Elaborando batch...');
tic;

batch_count = 0;
for i = 1:n_param1
    for j = 1:n_param2
        current_param = base_param;
        current_param.(param1_name) = param1_values(i);
        current_param.(param2_name) = param2_values(j);
        
        C1_result = popFlowV1MT_optimized(II, current_param);
        
        results{i, j, 1} = C1_result{1};
        results{i, j, 2} = C1_result{2};
        
        batch_count = batch_count + 1;
        
        if mod(batch_count, batch_size) == 0 || batch_count == total_combinations
            progress = batch_count / total_combinations;
            waitbar(progress, h, sprintf('Batch %d/%d completato', ...
                ceil(batch_count/batch_size), ceil(total_combinations/batch_size)));
            
            % Forza garbage collection ogni batch
            if mod(batch_count, batch_size*2) == 0
                clear ans;
                pause(0.1);
            end
        end
    end
end

close(h);
elapsed_time = toc;
fprintf('Completato in %.2f secondi\n', elapsed_time);

end