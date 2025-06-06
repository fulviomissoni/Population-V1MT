function memory_gb = estimate_memory_usage(sy, sx, n_frames, param, n_combinations)
% Stima utilizzo memoria approssimativo

n_orient = param.n_orient;
n_vel = length(param.pref_vel);

% Dimensioni matrice 5D base
base_elements = sy * sx * n_frames * n_orient * n_vel;

% Stima memoria per singola simulazione (in double precision)
single_sim_gb = base_elements * 8 * 10 / (1024^3); % fattore 10 per intermedi

% Memoria totale stimata
memory_gb = single_sim_gb * n_combinations * 0.5; % fattore sicurezza

end