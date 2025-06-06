% FUNZIONE DI ANALISI RISULTATI
function analyze_results(results, param1_name, param1_values, param2_name, param2_values)
% Analizza e visualizza i risultati del parameter sweep

fprintf('\n=== ANALISI RISULTATI ===\n');
fprintf('Parametro 1: %s (range: %.3f - %.3f)\n', param1_name, min(param1_values), max(param1_values));
fprintf('Parametro 2: %s (range: %.3f - %.3f)\n', param2_name, min(param2_values), max(param2_values));

[n_param1, n_param2, ~] = size(results);

% Calcola metriche per ogni combinazione
metrics = zeros(n_param1, n_param2, 4); % max, mean, std, energia_totale

for i = 1:n_param1
    for j = 1:n_param2
        for k = 1:2 % per entrambi gli output C1{1} e C1{2}
            data = results{i, j, k};
            if ~isempty(data)
                metrics(i, j, (k-1)*2 + 1) = max(data(:));
                metrics(i, j, (k-1)*2 + 2) = mean(data(:));
            end
        end
    end
end

% Trova combinazione ottimale (massima attivazione media)
mean_activation = mean(metrics(:,:,[2,4]), 3);
[max_val, linear_idx] = max(mean_activation(:));
[opt_i, opt_j] = ind2sub(size(mean_activation), linear_idx);

fprintf('\nCombinazione ottimale:\n');
fprintf('%s = %.4f, %s = %.4f\n', param1_name, param1_values(opt_i), param2_name, param2_values(opt_j));
fprintf('Attivazione media massima: %.6f\n', max_val);

% Visualizzazione opzionale
figure('Name', 'Parameter Sweep Results');
subplot(2,2,1);
imagesc(param2_values, param1_values, metrics(:,:,1));
title('Max Activation C1{1}');
xlabel(param2_name); ylabel(param1_name);
colorbar;

subplot(2,2,2);
imagesc(param2_values, param1_values, metrics(:,:,2));
title('Mean Activation C1{1}');
xlabel(param2_name); ylabel(param1_name);
colorbar;

subplot(2,2,3);
imagesc(param2_values, param1_values, metrics(:,:,3));
title('Max Activation C1{2}');
xlabel(param2_name); ylabel(param1_name);
colorbar;

subplot(2,2,4);
imagesc(param2_values, param1_values, metrics(:,:,4));
title('Mean Activation C1{2}');
xlabel(param2_name); ylabel(param1_name);
colorbar;

end