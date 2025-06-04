% make figures
close all
clc
clear all

load('SIMULATIONS/PlaidAnalysis/pop_resp_a2_0_5.mat')
addpath FUNCTIONS


sigma_r = 0.1;  
sigma_t = 0.4;  
K = 1.5; %inihibitory field size factor
% Thresholding parameter
 
logistic_slope = 8;
logistic_centre = 0.65;
max_iteration = 105;
%
theta_cell_OUT = 0:pi/param.n_orient:pi-pi/param.n_orient;    
[xx,tt] = meshgrid(param.pref_vel,theta_cell_OUT);
sigma = 0.25;
W2 = exp(-(xx(:).*cos(tt(:)-tt(:)') - xx(:)').^2./(2*sigma.^2)); %Andre pesi originali (inviluppo exp(.))
% W2 = exp(-(xx).^2/(2*sigma^2)).*cos(tt(:)-tt(:)'); %Chessa-Solari model
% W2 = cos(xx(:).*cos(tt(:)'-tt(:)) - xx(:)'); %primo modello per introdurre l'opponenza (inviluppo coseno)
% W2 = exp(-(xx(:).*cos(tt(:)'-tt(:)) - xx(:)').^2./(2*sigmaGabor.^2)).* ...
% cos(2*pi*1./(4*xx(:)').*(xx(:).*cos(tt(:)'-tt(:)) - xx(:)')); %inviluppo Gabor per gestire la quantità e la posizione dei pesi negativi
% W2 = reshape(reshape(W2,8,11,8,11)./max(reshape(W2,8,11,8,11),[],4),88,88);
W2(isnan(W2)) = 0;
pop_resp = reshape(pop_resp,param.n_orient,numel(param.pref_vel),numel(param.num_or_ch_pooled),numel(param.diff_c));
pop_resp(pop_resp<0) = 0;
PR_norm = pop_resp./mean(mean(pop_resp));
sz = size(PR_norm);
MT_norm = reshape(W2*reshape(PR_norm,sz(1)*sz(2),[]),sz);

stimIdx = param.diff_c == 0; % I will take a plaid made by gratings at same contrast levels
normIdx = 2;
[PR_decoded] = DecodeMxHat(squeeze(MT_norm(:,:,normIdx,stimIdx)),param,sigma_r,sigma_t,K,max_iteration,logistic_slope,logistic_centre);
fig=figure; plotPopResponse(squeeze(PR_decoded(:,:,max_iteration)),param.pref_vel)
set(gca,'Position',[0 0 1 1])
axis equal
colormap bone
title('Test Plaid - ')
% % saveas(gcf,"MT_P_popResp_PlaidII_Test_noNorm.pdf")
% saveas(gcf,"MT_P_popResp_PlaidII_Test_noNorm.png")
% savefig("MT_P_popResp_PlaidII_Test_noNorm")
% colorbar

stimIdx = param.diff_c == 0.5;
[PR_decoded] = DecodeMxHat(squeeze(MT_norm(:,:,normIdx,stimIdx)),param,sigma_r,sigma_t,K,max_iteration,logistic_slope,logistic_centre);
fig=figure; plotPopResponse(squeeze(PR_decoded(:,:,max_iteration)),param.pref_vel)
set(gca,'Position',[0 0 1 1])
colormap bone
axis equal
% % saveas(gcf,"MT_P_popResp_PlaidII_Ref_noNorm.pdf")
% saveas(gcf,"MT_P_popResp_PlaidII_Ref_noNorm.png")
% savefig("MT_P_popResp_PlaidII_Ref_noNorm")

% colorbar
%% psychometric
allcon = linspace(param.diff_c(1),param.diff_c(end-2),1000); 
y = 0.5 + (1-0.5)*(1-exp(-(allcon/0.5).^3));
figure, plot(flip(allcon),y,'k','LineWidth',2)
xlabel('\Deltac')
text(-0.13,1.04,'Pr(T|\theta,\Deltac_T,\Deltac_R)','FontSize',15)
set(gca,'Fontsize',15)

hold on
plot(allcon,0.75*ones(1,length(allcon)),'k--')
yticks([0.5,0.75,1])
% axis equal
set(gca,'Color','none')
box off
% 
% saveas(gcf,"psychometric_function",'fig')
% %export without background color
% fig2svg("pscychometric_function.svg",gcf)
% % saveas(gcf,"psychometric_function",'png')