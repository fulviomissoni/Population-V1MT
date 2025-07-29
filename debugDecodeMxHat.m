% Diagnostic version - let's see what happens at each step
function debugDecodeMxHat(pop_resp, param, sigma_r, sigma_t, K, max_iteration, logistic_slope, logistic_centre)

figure('Position', [100 100 1200 800]);

% Show original
subplot(2,4,1); 
imagesc(pop_resp); colorbar; 
title('Original Response');

% Build the kernel and show it
theta_cell_OUT = 0:pi/param.n_orient:pi-pi/param.n_orient;
[xx,tt] = meshgrid(param.pref_vel,theta_cell_OUT);
tt = tt + pi*(xx<0);
xx = abs(xx);
X = ((xx(:) - xx(:)').^2) / sigma_r^2 + ((tt(:) - tt(:)').^2) / sigma_t^2;
MX = 1/(2*pi*sigma_r*sigma_t)*exp(-X/2) - 1/(2*pi*K^2*sigma_r*sigma_t)*exp(-X/(2*K^2));
MX = MX./max(MX,[],'all');

subplot(2,4,2);
imagesc(reshape(MX(1,:), size(pop_resp))); colorbar;
title('Kernel (row 1)');

% Show first few iterations
th = 2e-2;
tmp = squeeze(pop_resp(:,:));
sze = size(tmp);

for iter = 1:min(4, max_iteration)
    CT = reshape(tmp,sze(1)*sze(2),1);
    CT(CT<th) = 0;
    M = max(CT,[],'all');
    
    % Show before logistic
    subplot(2,4,2+iter);
    imagesc(reshape(CT, sze)); colorbar;
    title(sprintf('Iter %d: Before Logistic (Max=%.3f)', iter, M));
    
    % Apply logistic
    if M > 0
        CT = M*1./(1+exp(-logistic_slope*(squeeze(CT)-M*logistic_centre)));
    end
    
    % Apply kernel
    CT = MX*CT;
    CT = CT./mean(CT,'all');
    
    tmp = reshape(CT, sze);
    
    if iter <= 2
        subplot(2,4,6+iter);
        imagesc(tmp); colorbar;
        title(sprintf('Iter %d: Final', iter));
    end
end
end