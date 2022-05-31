%
close all
clear 

addpath FUNCTIONS
load('SIMULATIONS\PlaidAnalysis\pop_resp.mat')
load('SIMULATIONS\MTWeights.mat')

my_pop_resp = pop_resp;
W3 = reshape(W,88,88)';
% W3 = W3./max(max(W3));
theta_cell_OUT = 0:pi/(param.n_orient):pi-pi/(param.n_orient);    
[xx,tt] = meshgrid(param.pref_vel,theta_cell_OUT);
% xx = [fliplr(xx(9:16,2:end)),xx(1:8,:)];
% tt = [fliplr(tt(9:16,2:end)),tt(1:8,:)];
sigma = 1;
W2 = exp(-(xx(:).*cos(tt(:)'-tt(:)) - xx(:)').^2/(2*sigma^2)).*cos(xx(:).*cos(tt(:)'-tt(:)) - xx(:)');
% W2 = exp(-((xx(:)'- xx(:)).^2)./(2*sigma^2)).*cos(tt(:)'-tt(:));
% W2 = W2 - eye(size(W2));
pop_resp_learning = reshape(W3*reshape(my_pop_resp,8*11,[]),8,11,2,11);
for i=1:11
    PR = squeeze(pop_resp_learning(:,:,1,i));
    PR(PR<0) = 0;
    figure, plotPopResponse(PR,param.pref_vel)
    colorbar
    pause
end