%MT PATTERN TUNING CURVES
close all
clear all
clc
addpath FUNCTIONS

load('SIMULATIONS/PlaidAnalysis/vel_tuning_All_PlaidII_lambda0.mat')
sze = size(EC);
idx_o = 1;
idx_v = 4;
%V1 tuning curves
PR_V1 = squeeze(EC(3,:,:,1,3,:,:));
TC_V1 = permute(PR_V1,[3 4 1 2]);
figure, plotPopResponse(squeeze(TC_V1(:,:,idx_o,idx_v)),param.prefVel)
title('V1 Tuning Curve')
figure, plotPopResponse(squeeze(PR_V1(:,:,idx_o,idx_v)),param.prefVel)
title('V1 Response')
%WEIGHTS
%Explicit intersection of constraints method to compute weigths
theta_cell_OUT = 0:pi/param.nOrient:pi-pi/param.nOrient;    
[xx,tt] = meshgrid(param.prefVel,theta_cell_OUT);
% sigmaGabor = 0.5*4/3*xx(:)'; %sigma = 0.3*4/3*xx(:);
sigma = 0.25;
W2 = exp(-(xx(:).*cos(tt(:)'-tt(:)) - xx(:)').^2./(2*sigma.^2));
% W2 = exp(-(xx(:)-xx(:)').^2./(2*sigma.^2)).*cos(tt(:)'-tt(:));
% W2 = cos(xx(:).*cos(tt(:)'-tt(:)) - xx(:)');
% W2 = exp(-(xx(:).*cos(tt(:)'-tt(:)) - xx(:)').^2./(2*sigmaGabor.^2)).*cos((xx(:).*cos(tt(:)'-tt(:)) - xx(:)'));
% W2 = exp(-(xx(:).*cos(tt(:)'-tt(:)) - xx(:)').^2./(2*sigmaGabor.^2)).*cos(2*pi*1./(4*xx(:)').*(xx(:).*cos(tt(:)'-tt(:)) - xx(:)'));
% W2 = reshape(reshape(W2,8,11,8,11)./max(reshape(W2,8,11,8,11),[],4),88,88);
W2(isnan(W2)) = 0;
%WEIGHTS
Wplot = reshape(W2,8,11,8,11);
figure, plotPopResponse(squeeze(Wplot(idx_o,idx_v,:,:)),param.prefVel)
title('WEIGHTS')

%MT tuning curves
PR_MT = reshape(W2*reshape(PR_V1,88,88),8,11,8,11);
TC_MT = permute(PR_MT,[3,4,1,2]);
figure, plotPopResponse(squeeze(TC_MT(:,:,idx_o,idx_v)),param.prefVel)
title('MT Pattern Tuning Curve')
figure, plotPopResponse(squeeze(PR_MT(:,:,idx_o,idx_v)),param.prefVel)
title('MT response')
%% reorganize to have 16 orientation channels

myTC_MT = squeeze(TC_MT(:,:,idx_o,idx_v));
myTC_MT16 = [myTC_MT(:,6:end);fliplr(myTC_MT(:,1:6))];
figure, plot(rad2deg(0:pi/param.nOrient:2*pi-pi/param.nOrient),squeeze(myTC_MT16(:,idx_v)))
title('MT Tuning curve along one velocity channel')
%Pool = [fliplr(PFull(9:16,2:end)),PFull(1:8,:)];