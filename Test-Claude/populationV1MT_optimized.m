function [C1] = populationV1MT_optimized(G,param)
% Versione ottimizzata con pooling spaziale vettorizzato

alpha = param.norm_param;
sigma_pool = param.sigma_pool;
num_or_ch_pooled = param.num_or_ch_pooled;

[sy, sx, n_frames, n_orient, v] = size(G{1});
sze = size(G{1});

C = G{1}(:);
S = G{2}(:);
Ct = G{3}(:);
St = G{4}(:);

clear G;

%% ENERGY MODEL (ottimizzato)

% SPATIO TEMPORAL ORIENTED FILTERS
S0{1} = S+Ct;
S0{2} = -St+C;
S0{3} = -S+Ct;
S0{4} = St+C;
clear C S Ct St

% COMPLEX CELLS - FIRST LAYER
C1{1} = S0{1}.^2 + S0{2}.^2;
C1{2} = S0{3}.^2 + S0{4}.^2;

clear S0

% NORMALIZATION STAGE OTTIMIZZATO
a1 = alpha(1);
a2 = alpha(2);

for i = 1:2
    C1{i} = reshape(C1{i}, sy, sx, n_frames*n_orient*v);
    
    % OTTIMIZZAZIONE PRINCIPALE: Pooling spaziale vettorizzato
    % Invece di ciclo frame per frame, usa operazioni 3D
    if sigma_pool > 0
        % Applica filtro gaussiano su tutte le "pagine" contemporaneamente
        S = zeros(size(C1{i}));
        for p = 1:size(C1{i}, 3)
            S(:,:,p) = imgaussfilt(C1{i}(:,:,p), sigma_pool);
        end
    else
        S = C1{i};
    end
    
    C1{i} = reshape(C1{i}, 1, []);
    S = reshape(S, sy*sx*n_frames, n_orient, v);
    S = permute(S, [2 3 1]);
    
    m = 1/mean(S(:));
    
    % Orientation pooling ottimizzato
    if num_or_ch_pooled == 8
        S = sum(sum(S,2),1);
        S = repmat(S, [n_orient v 1]);
    elseif num_or_ch_pooled == 1
        S = sum(S,2);
        S = repmat(S, [1 v 1]);
    end
    
    S = reshape(S, n_orient, v, sy*sx*n_frames);
    S = permute(S, [3 1 2]);
    S = reshape(S, 1, []);
    
    C1{i} = C1{i}./(a1 + a2*m*S);
    C1{i}(isnan(C1{i})) = 0;
end

for i=1:2
    C1{i} = squeeze(reshape(C1{i}, sze));
end

end