function C1 = popFlowV1MT_optimized(II, param)
% Versione ottimizzata della funzione originale
% Ottimizzazioni principali:
% 1. Vettorizzazione del pooling spaziale
% 2. Pre-allocazione ottimizzata
% 3. Riduzione copie di memoria

k0 = param.spat_freq;
v = param.pref_vel;
Ft_choice = param.temp_filt;
filter_file = param.spatial_filt;

% COMPUTE GABOR FILTERING (mantenuto originale)
if param.n_orient == 8
    [F] = filtGaborSpace2D(II,filter_file);
end
if param.n_orient == 16
    [F] = filtSepGabor(II,filter_file);
end

F = filtTime(F,'valid',Ft_choice,v,k0);

for i=1:4
    F{i} = permute(F{i},[1 2 4 3 5]);
end

% POPULATION ENCODING OTTIMIZZATO
C1 = populationV1MT_optimized(F,param);

clear F

end