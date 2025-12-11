function W = weightingfunctions_general(xx, tt)
% Generalized weighting function for orientation channels in [0, pi]
% with implicit 2*pi periodicity handled via negative velocities
%
% Inputs:
%   xx - velocity values (column vector, length = n_orient * n_vel)
%   tt - orientation values (column vector, length = n_orient * n_vel)
%   n_orient - number of orientation channels (e.g., 8)
%   n_vel - number of velocity channels (e.g., 11)
%
% Outputs:
%   W - weighting matrices [N x N x 5] where N = n_orient * n_vel

sigma_t = 1;
N = numel(xx);
W = zeros([N, N, 5]);

% Weight 1: Gaussian envelope
W(:,:,1) = exp(-(xx(:).*cos(tt(:)-tt(:)') - xx(:)').^2 ./ (2*sigma_t^2));

tt_eff = tt + pi * (xx < 0);
% Weight 2: Cosine tuning with 2*pi periodicity
% cos(theta_preferred_output - theta_input)
W(:,:,2) = cos(tt_eff(:) - tt_eff(:)');

% Weight 3: Cosine opponent model
W(:,:,3) = cos(xx(:).*cos(tt(:)'-tt(:)) - xx(:)');

% Weight 4: Gabor envelope
gabor_arg = xx(:).*cos(tt(:)'-tt(:)) - xx(:)';
W(:,:,4) = exp(-gabor_arg.^2 ./ (2*sigma_t^2)) .* cos(2*pi*gabor_arg ./ (4*xx(:)'));

% Weight 5: Andre original weights (linear)
W(:,:,5) = xx(:).*cos(tt(:)-tt(:)') - xx(:)';

W(isnan(W)) = 0;

end