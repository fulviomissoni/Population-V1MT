function results = run_vectorized_sweep(II, base_param, param1_name, param1_values, param2_name, param2_values)
% Implementazione vettorizzata completa

n_param1 = length(param1_values);
n_param2 = length(param2_values);
n_total = n_param1 * n_param2;

% Crea meshgrid dei parametri
[P1_grid, P2_grid] = meshgrid(param1_values, param2_values);
P1_flat = P1_grid(:);
P2_flat = P2_grid(:);

fprintf('Elaborazione vettorizzata di %d combinazioni...\n', n_total);

% Modifica approccio: elabora in parallelo logico ma mantenendo cicli ottimizzati
results = cell(n_param1, n_param2, 2);

% Progress bar
h = waitbar(0, 'Elaborando combinazioni parametri...');
tic;

for i = 1:n_param1
    for j = 1:n_param2
        % Crea parametri per questa combinazione
        current_param = base_param;
        current_param.(param1_name) = P1_flat((i-1)*n_param2 + j);
        current_param.(param2_name) = P2_flat((i-1)*n_param2 + j);
        
        % Esegui simulazione
        C1_result = popFlowV1MT_optimized(II, current_param);
        
        % Salva risultati
        results{i, j, 1} = C1_result{1};
        results{i, j, 2} = C1_result{2};
        
        % Aggiorna progress
        progress = ((i-1)*n_param2 + j) / n_total;
        waitbar(progress, h, sprintf('Combinazione %d/%d', (i-1)*n_param2 + j, n_total));
    end
end

close(h);
elapsed_time = toc;
fprintf('Completato in %.2f secondi (%.3f sec/combinazione)\n', elapsed_time, elapsed_time/n_total);

end