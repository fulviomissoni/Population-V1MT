function results = run_parameter_sweep_vectorized(II, base_param, param1_name, param1_values, param2_name, param2_values)
% Esegue parameter sweep vettorizzato per popFlowV1MT
% 
% Input:
%   II - immagini input [rows, cols, frames]
%   base_param - struttura parametri base
%   param1_name, param2_name - nomi dei parametri da variare (string)
%   param1_values, param2_values - vettori dei valori da testare
%
% Output:
%   results - cell array 4D: {param1_idx, param2_idx, output_idx}

fprintf('Iniziando parameter sweep: %d x %d = %d combinazioni\n', ...
    length(param1_values), length(param2_values), ...
    length(param1_values) * length(param2_values));

% Pre-allocazione risultati
n_param1 = length(param1_values);
n_param2 = length(param2_values);
results = cell(n_param1, n_param2, 2); % 2 per i due output C1{1} e C1{2}

% Stima memoria necessaria
[sy, sx, n_frames] = size(II);
estimated_memory_gb = estimate_memory_usage(sy, sx, n_frames, base_param, n_param1 * n_param2);
fprintf('Memoria stimata necessaria: %.2f GB\n', estimated_memory_gb);

if estimated_memory_gb > 16
    fprintf('ATTENZIONE: Memoria elevata. Utilizzando elaborazione a batch.\n');
    results = run_batched_sweep(II, base_param, param1_name, param1_values, param2_name, param2_values);
    return;
end

% APPROCCIO 1: Vettorizzazione completa (se memoria sufficiente)
try
    results = run_vectorized_sweep(II, base_param, param1_name, param1_values, param2_name, param2_values);
catch ME
    if contains(ME.message, 'memory') || contains(ME.message, 'Out of memory')
        fprintf('Memoria insufficiente per vettorizzazione completa. Utilizzando batch processing.\n');
        results = run_batched_sweep(II, base_param, param1_name, param1_values, param2_name, param2_values);
    else
        rethrow(ME);
    end
end

end