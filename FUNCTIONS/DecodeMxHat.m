function [pop_resp_V1MT, v_CM, v_max] = DecodeMxHat( ...
    pop_resp, param, sigma_r, sigma_t, K, max_iteration, logistic_slope, logistic_centre)
%DECODEMXHAT Competitive localization via Mexican-hat filtering in neural space

%% normalise pop_resp for having weights for CM decoding
% Normalize pop_resp to ensure weights for CM decoding
pop_resp = pop_resp / sum(pop_resp(:));
%% Parameters
% th = 2e-2;
th = 1e-15;
theta_cell_OUT = 0:pi/param.n_orient:pi-pi/param.n_orient;

n_orient = param.n_orient;
n_vel = numel(param.pref_vel);

%% ================= INTERPOLATION (unchanged) =================
[tt_nd, xx_nd] = ndgrid(theta_cell_OUT, param.pref_vel);

tt_fine = linspace(min(theta_cell_OUT), max(theta_cell_OUT), 100);
xx_fine = linspace(min(param.pref_vel), max(param.pref_vel), 100);
[tt_interp, xx_interp] = ndgrid(tt_fine, xx_fine);

pop_resp_interp = interp2(xx_nd, tt_nd, pop_resp, xx_interp, tt_interp, 'cubic');

% tt_interp_adj = tt_interp + pi*(xx_interp < 0);
% xx_interp_adj = abs(xx_interp);



%% ================= MEXICAN HAT (FIXED) =================
% Neural grid (orientation × velocity)
[TH1, V1] = ndgrid(theta_cell_OUT, param.pref_vel);
[TH2, V2] = ndgrid(theta_cell_OUT, param.pref_vel);

TH1 = TH1(:);  V1 = V1(:);
TH2 = TH2(:);  V2 = V2(:);

% --- circular orientation distance (π-periodic) ---
dtheta = abs(TH1 - TH2');
dtheta = min(dtheta, pi - dtheta);

% --- velocity distance ---
dvel = abs(V1 - V2');

% --- squared distance in NEURAL space ---
D2 = (dvel.^2) / sigma_r^2 + (dtheta.^2) / sigma_t^2;

% --- Mexican hat (Difference of Gaussians) ---
Exc = exp(-D2/2);
Inh = (1/K^2) * exp(-D2/(2*K^2));
MX = Exc - Inh;

% --- enforce competition (zero-sum rows) ---
% MX = MX - mean(MX, 2);

%% ================= ITERATIVE DYNAMICS (CORRECTED) =================
pop_resp_V1MT = zeros(n_orient, n_vel, max_iteration);
pop_resp_V1MT(:,:,1) = pop_resp;

tmp = pop_resp(:);

for it = 2:max_iteration
    % 1. Thresholding (Noise suppression)
    tmp(tmp < th) = 0;
    
    % 2. Apply Filter
    drive = MX * tmp;
    
    % 3. Apply Update (Rectified)
    % We add the drive to the previous state (integration) or replace it.
    % For simple filtering, replacement is fine:
    tmp = drive;
    tmp(tmp < 0) = 0; % Ensure non-negative after inhibition
    
    % % 4. NORMALIZATION (Crucial Replacement for Logistic)
    % % This prevents explosion and allows the "blob" to stabilize.
    % % Dividing by the max value keeps the peak height = 1, 
    % % allowing the width to settle based on sigma_r/sigma_t.
    % current_max = max(tmp(:));
    % if current_max > 0
    %     tmp = tmp ./ current_max; 
    % end
    
    % Store
    pop_resp_V1MT(:,:,it) = reshape(tmp, n_orient, n_vel);
end

%% ================= CENTER OF MASS DECODING (INTERPOLATED) =================

% Risposta finale
subPopResp = pop_resp_V1MT(:,:,end)./ sum(pop_resp_V1MT(:,:,end),'all');

% --- interpolazione della risposta finale ---
subPopResp_interp = interp2( ...
    xx_nd, tt_nd, subPopResp, ...
    xx_interp, tt_interp, ...
    'cubic', 0);

% % --- stessa correzione angolo/velocità usata per v_max ---
% tt_interp_adj = tt_interp + pi*(xx_interp < 0);
% xx_interp_adj = abs(xx_interp);

% --- centro di massa nello spazio interpolato ---
% M = sum(subPopResp_interp(:));

v_CM(1) = sum(sum(subPopResp_interp .* ...
           (xx_interp .* cos(tt_interp)))) ./ sum(sum(subPopResp_interp));

v_CM(2) = sum(sum(subPopResp_interp .* ...
           (xx_interp .* sin(tt_interp)))) ./ sum(sum(subPopResp_interp));

[~, max_ind] = max(subPopResp_interp(:));
[row_max, col_max] = ind2sub(size(subPopResp_interp), max_ind);

v_max(1) = xx_interp(row_max,col_max) * cos(tt_interp(row_max,col_max));
v_max(2) = xx_interp(row_max,col_max) * sin(tt_interp(row_max,col_max));

end
