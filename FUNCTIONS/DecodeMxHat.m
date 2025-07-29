function [pop_resp_V1MT, vx, vy] = DecodeMxHat(pop_resp,param,sigma_r,sigma_t,K,max_iteration,logistic_slope,logistic_centre)
%DECODEMXHAT Iterative competitive localization of neural population activity
%
%   [pop_resp_V1MT, vx, vy] = DECODEMXHAT(pop_resp, param, sigma_r, sigma_t, K, max_iteration, logistic_slope, logistic_centre)
%   performs iterative competitive localization on neural population responses using
%   center-surround filtering with recurrent processing and nonlinear transformations.
%   The algorithm enhances local activity peaks while suppressing diffuse background
%   responses, then decodes the resulting pattern using population vector decoding.
%
%   INPUTS:
%   pop_resp        - Population response matrix [n_orient x n_vel]
%                     Neural activity where rows=orientation channels, cols=velocity preferences
%   param           - Parameter structure with required fields:
%                     .n_orient: number of orientation channels (integer)
%                     .pref_vel: preferred velocity values (vector)  
%   sigma_r         - Velocity dimension spatial scale (positive scalar, typically 0.1-10)
%                     Controls extent of center-surround filtering in velocity space
%   sigma_t         - Orientation dimension spatial scale (positive scalar, typically 0.1-π/2)  
%                     Controls extent of center-surround filtering in orientation space
%   K               - Inhibition scale factor (scalar > 1, typically 1.5-5)
%                     Ratio of inhibitory to excitatory spatial spread
%                     K=1.5: strong lateral inhibition; K>3: gentle inhibition
%   max_iteration   - Maximum number of recurrent iterations (integer, typically 5-15)
%                     More iterations = stronger localization and noise suppression
%   logistic_slope  - Logistic function steepness (positive scalar, typically 1-20)
%                     Controls aggressiveness of contrast enhancement  
%   logistic_centre - Logistic function threshold (scalar 0-1, typically 0.1-0.5)
%                     Threshold level as fraction of maximum response
%
%   OUTPUTS:
%   pop_resp_V1MT   - Processed population response [n_orient x n_vel x max_iteration]
%                     Activity after each iteration; use (:,:,end) for final result
%   vx              - Decoded x-component of velocity (scalar, same units as pref_vel)
%                     Center-of-mass velocity from final processed response
%   vy              - Decoded y-component of velocity (scalar, same units as pref_vel)
%                     Center-of-mass velocity from final processed response
%
%   ALGORITHM:
%   The function implements a 4-stage iterative process:
%   1. Distance Computation: Creates 2D feature space combining velocity and orientation
%   2. Center-Surround Kernel: Builds difference-of-Gaussians with gentle inhibition
%   3. Iterative Processing: Applies recurrent competitive dynamics:
%      - Hard thresholding (noise removal)
%      - Logistic contrast enhancement
%      - Center-surround filtering  
%      - Mean normalization (stability)
%   4. Population Vector Decoding: Center-of-mass from final activity pattern
%
%   TYPICAL USAGE:
%   % Standard motion processing parameters
%   [processed, vx, vy] = DecodeMxHat(pop_resp, param, 2.0, pi/8, 2.5, 10, 5, 0.2);
%
%   % For fine localization (stronger competition)
%   [processed, vx, vy] = DecodeMxHat(pop_resp, param, 1.0, pi/8, 1.8, 15, 10, 0.2);
%
%   % For noise robustness (gentler processing)  
%   [processed, vx, vy] = DecodeMxHat(pop_resp, param, 3.0, pi/6, 3.5, 8, 3, 0.3);
%
%   NOTES:
%   - Algorithm preserves multiple activity peaks unlike winner-take-all approaches
%   - Biologically plausible: mimics recurrent cortical center-surround dynamics
%   - Robust to noise through iterative refinement
%   - Competition strength tunable via K parameter
%   - Convergence typically achieved within 10-15 iterations
%
%   See also: Population vector decoding, center-surround filtering, competitive networks
%
%   Author: Fulvio Missoni; Andrea Canessa;
%   Date: 28-07-25
%   Version: 1.0


%%
% Set parameters
th = 2e-2;
diff_c = param.diff_c;
num_or_ch_pooled = param.num_or_ch_pooled;

theta_cell_OUT = 0:pi/param.n_orient:pi-pi/param.n_orient;
[xx,tt] = meshgrid(param.pref_vel,theta_cell_OUT);
tt = tt + pi*(xx<0);
xx = abs(xx);
dx = xx(:).*cos(tt(:));
dy = xx(:).*sin(tt(:));
X = ((xx(:) - xx(:)').^2) / sigma_r^2 + ((tt(:) - tt(:)').^2) / sigma_t^2;
%MEXICAN HAT
MX = 1/(2*pi*sigma_r*sigma_t)*exp(-X/2) - ...
    1/(2*pi*K^2*sigma_r*sigma_t)*exp(-X/(2*K^2));
MX = MX./max(MX,[],'all');

%Max number of recurrency iteration
%     max_iteration = 10;
sze = size(squeeze(pop_resp(:,:)));

% logistic_centre = M*logistic_centre;
% for i=1:size(pop_resp,1)
    % tmp = pop_resp_BioGautama;
    tmp = squeeze(pop_resp(:,:));
    %organize population responses
    pop_resp_V1MT(:,:,1) = squeeze(pop_resp(:,:));
    for indResp = 2:max_iteration
        %Apply my Weights
        %iterates mexican hat weigthing function
        CT = reshape(tmp,sze(1)*sze(2),1);
        CT(CT<th) = 0;
        M = max(CT,[],'all');
        CT = M*1./(1+exp(-logistic_slope*(squeeze(CT)-M*logistic_centre)));
        CT = MX*CT;
        CT = CT./mean(CT,'all');
        pop_resp_V1MT(:,:,indResp) = reshape(CT,sze);
        tmp = CT;
    end
% end
%Thresholding
% M = max(squeeze(pop_resp_V1MT(:,:,max_iteration)),[],'all');
pop_resp_V1MT(:,:,max_iteration) = M*1./(1+exp(-logistic_slope*(squeeze(pop_resp_V1MT(:,:,max_iteration))-M*logistic_centre)));
pop_resp_V1MT = reshape(pop_resp_V1MT,param.n_orient,numel(param.pref_vel),max_iteration);

% centre of mass
% for i=1:numel(diff_c)
subPopResp = squeeze(pop_resp_V1MT(:,:,max_iteration));
M = sum(sum(subPopResp));
%centre of mass decoding
vx = squeeze(sum(sum(subPopResp.*(xx(1:param.n_orient,:).*cos(tt(1:param.n_orient,:))),1),2)./M);
vy = squeeze(sum(sum(subPopResp.*(xx(1:param.n_orient,:).*sin(tt(1:param.n_orient,:))),1),2)./M);
end
