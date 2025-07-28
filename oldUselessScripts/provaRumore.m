clear all
% close all
clc
addpath FUNCTIONS

load('SIMULATIONS/PlaidAnalysis/pop_resp_a2_0_5.mat')
load('SIMULATIONS/MTWeights.mat')
theta_cell_OUT = 0:pi/param.n_orient:pi-pi/param.n_orient;    
[xx,tt] = meshgrid(param.pref_vel,theta_cell_OUT);

sigma_r = 0.1;  
sigma_t = 0.4;  
K = 1.5; %inihibitory field size factor
% Thresholding parameter
 
logistic_slope = 1.5;
logistic_centre = 0.6;
max_iteration = 3;

%     W = W/max(W(:));
%Explicit intersection of constraints method to compute weigths
% sigmaGabor = 0.5*4/3*xx(:)'; %sigma = 4/3*xx(:)';
sigma = 0.25;
W2 = exp(-(xx(:).*cos(tt(:)-tt(:)') - xx(:)').^2./(2*sigma.^2)); %Andre pesi originali (inviluppo exp(.))
% W2 = exp(-(xx).^2/(2*sigma^2)).*cos(tt(:)-tt(:)'); %Chessa-Solari model
% W2 = cos(xx(:).*cos(tt(:)'-tt(:)) - xx(:)'); %primo modello per introdurre l'opponenza (inviluppo coseno)
% W2 = exp(-(xx(:).*cos(tt(:)'-tt(:)) - xx(:)').^2./(2*sigmaGabor.^2)).* ...
% cos(2*pi*1./(4*xx(:)').*(xx(:).*cos(tt(:)'-tt(:)) - xx(:)')); %inviluppo Gabor per gestire la quantità e la posizione dei pesi negativi
% W2 = reshape(reshape(W2,8,11,8,11)./max(reshape(W2,8,11,8,11),[],4),88,88);
W2(isnan(W2)) = 0;

% W3 = reshape(W,88,88);
pop_resp = reshape(pop_resp,param.n_orient,numel(param.pref_vel),numel(param.num_or_ch_pooled),numel(param.diff_c));
pop_resp(pop_resp<0) = 0;
PR_norm = pop_resp./mean(mean(pop_resp));
sz = size(PR_norm);
MT_norm = reshape(W2*reshape(PR_norm,sz(1)*sz(2),[]),sz);

% logistic_centre = 0.7*M;
% PR = reshape(PR,param.n_orient,numel(param.pref_vel),numel(param.num_or_ch_pooled),numel(param.diff_c));

% pop_resp = reshape(pop_resp,param.n_orient,numel(param.pref_vel),numel(param.num_or_ch_pooled),numel(param.diff_c));
% sze = size(pop_resp);
% PR_noNorm = reshape(W3*reshape(squeeze(pop_resp(2,1,:,:)),param.n_orient*numel(param.pref_vel),[]),param.n_orient,numel(param.pref_vel));
% PR_norm = reshape(W3*reshape(squeeze(pop_resp(1,1,:,:)),param.n_orient*numel(param.pref_vel),[]),param.n_orient,numel(param.pref_vel));%WEIGHTS

% pop_resp_BioGautama = reshape(pop_resp_BioGautama,param.n_orient,numel(param.pref_vel),numel(param.num_or_ch_pooled),numel(param.diff_c));
%%
mySeed = 1;
rng(mySeed)
nTrial = 1000;
M = mean(mean(mean(mean(MT_norm))));

refIdx = find(param.diff_c==0.8);
sigma_ref = 2e-1*mean(std(mean(MT_norm(:,:,:,refIdx),3)));

for n=1:nTrial
for j=1:numel(param.num_or_ch_pooled)
    for i=1:numel(param.diff_c)
        Test_tmp = squeeze(MT_norm(:,:,j,i));
        Ref_tmp = squeeze(MT_norm(:,:,j,refIdx));        
        %THRESHOLDING
%         Test_tmp = M*logistic(squeeze(MT_norm(:,:,j,i)),logistic_centre,logistic_slope);
%         Ref_tmp = M*logistic(squeeze(MT_norm(:,:,j,9)),logistic_centre,logistic_slope);
%         Test_tmp = reshape(W2*reshape(squeeze(pop_resp(:,:,j,i)),88,[]),8,11);
%         Ref_tmp = reshape(W2*reshape(squeeze(pop_resp(:,:,j,9)),88,[]),8,11);  
    
        sigma = 2e-1*mean(std(Test_tmp)); 
        Test_tmp = Test_tmp + sigma*randn(size(Test_tmp));
%         Test_tmp = Test_tmp + poissrnd(mean(mean(Test_tmp)));
        PR.T(:,:,j,i) = Test_tmp;
%         PR.T(:,:,j,i) = reshape(W2*reshape(tmp,param.n_orient*numel(param.pref_vel),[]),param.n_orient,numel(param.pref_vel));

%         sigma = 2e-1*mean(std(Ref_tmp)); 
        Ref_tmp = Ref_tmp + sigma_ref*randn(size(Ref_tmp));
%         Ref_tmp = Ref_tmp + poissrnd(mean(mean(Ref_tmp)));
        PR.R(:,:,j,i) = Ref_tmp;
%         PR.R(:,:,j,i) = reshape(W2*reshape(tmp,param.n_orient*numel(param.pref_vel),[]),param.n_orient,numel(param.pref_vel));

        %DECODING ACTIVITY, FIRST METHOD - take average activity along truetheta direction
        %Select truetheta and evaluate pop_activity along that direction
        PR_truetheta(:,1) = squeeze(Ref_tmp(3,:));
        PR_truetheta(:,2) = squeeze(Test_tmp(3,:));
        
        values(:,i,n) = mean(PR_truetheta);
        [dump, resp(i)] = max(mean(PR_truetheta));
        %DECODING ACTIVITY, SECOND METHOD - CM DECODING
        [PR.R_noise(:,:,:,i,j), vx(i,j,n,1),vy(i,j,n,1)] = DecodeMxHat(squeeze(Ref_tmp),param,sigma_r,sigma_t,K,max_iteration,logistic_slope,logistic_centre);
        [PR.T_noise(:,:,:,i,j), vx(i,j,n,2),vy(i,j,n,2)] = DecodeMxHat(squeeze(Test_tmp),param,sigma_r,sigma_t,K,max_iteration,logistic_slope,logistic_centre);
    end
%     figure, plot(param.diff_c,rad2deg(stim.truetheta-squeeze(atan2(vy(:,j,:),vx(:,j,:)))))
%     legend('Reference','Test')
%     xlabel('\Deltac')
%     title(['numOrPooled =',num2str(param.num_or_ch_pooled(j))])
%     ylabel('\Delta\theta [degrees]')
%     ylim([0,60])
%     figure, plot(param.diff_c,diff(values))
%     hold on, plot(param.diff_c,zeros(1,numel(param.diff_c)),'r--')
%     legend('En_{Test}-En_{Ref}','Threshold')
%     xlabel('\Deltac')
%     title(['numOrPooled =',num2str(param.num_or_ch_pooled(j))])
%     figure, bar(1:2,values')
%     ylabel('mean activity')
%     xticks([1,2])
%     xticklabels({'Reference','Test'})
%     title(['numOrPooled =',num2str(param.num_or_ch_pooled(j))])
%     legend(num2str(param.diff_c'))
end
end
%%
%Normalization effect on pop activity on unbalanced plaid perception
figure, plotPopResponse(squeeze(PR.R_noise(:,:,max_iteration,5,1)),param.pref_vel)
colorbar
title('Reference No Normalized')
figure, plotPopResponse(squeeze(PR.T_noise(:,:,max_iteration,5,1)),param.pref_vel)
colorbar
title('Test No Normalized')
figure, plotPopResponse(squeeze(PR.R_noise(:,:,max_iteration,refIdx+1,2)),param.pref_vel)
title('Reference Normalized')
figure, plotPopResponse(squeeze(PR.T_noise(:,:,max_iteration,refIdx+1,2)),param.pref_vel)
title('Test Normalized')
% %% what do you choose?
% angle_err = stim.truetheta-atan2(vy,vx);
% risp = permute(angle_err,[1 3 2]);
% risp = squeeze(abs((risp(:,2,:)-risp(:,1,:))).^2<1e-3);
% figure, plot(param.diff_c,risp)

%% mean and standard deviation plot
theta_est = stim.truetheta - atan2(vy,vx);
%Eight Or Pooled ==1 ||One Or Pooled ==2
for normIdx = 1:2
    stimIdx = 1; %Ref==1, Test==1
    
    figure, PlotMeanStd(param.diff_c,rad2deg(squeeze(theta_est(:,normIdx,:,stimIdx)))','b')
    stimIdx = 2; %Ref==1, Test==1
  
    hold on, PlotMeanStd(param.diff_c,rad2deg(squeeze(theta_est(:,normIdx,:,stimIdx)))','r')
    ylim([0,50])
end
%% Psychometric function
myText = ["no",""];
xlabels = ["0","0.2","0.4","0.6","0.8"];
filename = strcat("psychometric_model");
x = flip(param.diff_c(1:9));
figure, 
for normIdx = 1:2
    risp = squeeze(theta_est(:,normIdx,:,1)-theta_est(:,normIdx,:,2)) > 0;
    pr = squeeze(sum(risp(3:end,:),2)/nTrial);
%     pr(pr<0.5)=0;
    sz = 80; % [];
    c = zeros(3,length(x))+0.7;
    indexes = pr>=0.5;
    tmp = pr(indexes);
    scatter(x(indexes),tmp,sz,c(:,indexes)','filled')
    ylim([0.5,1])

    colormap gray

%     figure, plot(param.diff_c,pr')
%     ft = fittype('0.5 + (1-0.5)*(1-exp(-(x/a).^b))'); % Weibull function for 2AFC
    ft = fittype( @(a,b,x)Weibull(a,b,x));
    x_tofit = x;
    y_tofit = pr;
    ft_fit = fit(x_tofit',y_tofit, ft,'Lower',[0,0],'Upper',[1,Inf],'StartPoint',[.1 1]);
    allcon = linspace(x(1),x(end),1000); 
    y_fit = 0.5 + (1-0.5)*(1-exp(-(allcon/ft_fit.a).^ft_fit.b));
                
    hold on
    if normIdx == 1
        cdata = [0.2,0.2,0.2];
    else
        cdata = [1,0,0];
    end
    plot(allcon,y_fit,'color',cdata,'linew',2)
    hold on
    %threshold
    idx = find(y_fit<0.75);
    idx = idx(1);
    plot(flip(allcon(idx:end)),0.75*ones(1,numel(allcon(idx:end))),'k--')
    plot(allcon(idx)*ones(1,numel(allcon(1:idx))),linspace(0.5,y_fit(idx),numel(allcon(1:idx))),'k--')

    scatter(allcon(idx),y_fit(idx),30,[0.2,0.2,0.2],'filled')
    ylabel('Pr(T|\theta,\Deltac_T,\Deltac_R)')
    xlabel('\Deltac')
%     
%     xlim([0,0.8])
    ylim([0.5,1])
    set(gca,'Fontsize',15)

    set(gca,'Color','none')
    box off

%     title(['Psychometric Function ',myText(normIdx),' norm'])
%     ylim([0,1.2])
%     figure, PlotMeanStd(param.diff_c,pr')
% %     ylim([0,1.2])
%     xlim([0,0.8])
end
% text(-0.08,1.025,'Pr(T|\theta,\Deltac_T,\Deltac_R)','FontSize',15)
% xlim([0,0.8])
xticks(flip(x(1:2:end))) 
xticklabels(flip(xlabels));
ylim([0.5,1])
fig2svg(strcat(filename,".svg"))
legend('','L=8','','','','','L=1','','','')
savefig(strcat(filename,".fig"))