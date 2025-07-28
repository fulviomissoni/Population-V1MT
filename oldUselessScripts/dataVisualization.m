clear 
close all
clc

addpath FUNCTIONS

%% IT DOESN'T WORK
% %% load data and organize them in cell array
% % lambda = ["0","1e-06","1e-05","00001","00001","0001","001"];
% lambda = ["0","1e-06","1e-05","1e-04","1e-03"];
% for eta=1:numel(lambda)
% % block = 'vel_tuning_All_PlaidII_lambda';
% block = 'newVel_tuning_PlaidII_lambda';
% mypath = 'SIMULATIONS\';
% filename = [mypath,block,convertStringsToChars(lambda(eta)),'_difContrasts.mat'];
% % filename = [mypath,block,convertStringsToChars(lambda(1)),'.mat'];
% load(filename)
% %%
% %size(e) is known
% nTheta_g = [6,7];
% nVel = 5; %vel stim tested -2pix/frames and -1.2 pix/frames
% nContr_g = 4; %contrast stim tested are grat1 = [0,0.2,0.4]; grat2 = [1,0.8,0.6];
% pop_resp{1} = squeeze(e(1,:,:,:));
% theta_cell_OUT = 0:pi/param.nOrient:pi-pi/param.nOrient;
% [xx,tt] = meshgrid(param.prefVel,theta_cell_OUT);
% sze = [param.nOrient,numel(param.prefVel),nTheta_g(1),nTheta_g(2),nVel,nContr_g];
% [rho, theta] = meshgrid(param.prefVel,linspace(0,pi,sze(1)+1));
% 
% % LEARNED WEIGHTS
% load([mypath,'GautamaWeights88_Plaid.mat'])
% % EXPLICIT IOC WEIGHTS
% W(:,:,2) = exp(-(xx(:).*cos(tt(:)'-tt(:)) - xx(:)').^2/(2*0.25^2)); % my weights
% W(:,:,2) = W(:,:,2)-W(:,:,2).*eye(size(W(:,:,2)));
% % MEXICAN HAT WEIGHTS
% %parameters
% sigma_r = 0.1;  
% sigma_t = 0.4;  
% K = 1.5; %inihibitory field size factor
% %Thresholding function parameter
% slopeLogistic = 9;  
% centerLogistic = 0.7;
% 
% [xx,tt] = meshgrid(param.prefVel,theta_cell_OUT);
% tt = tt + pi*(xx<0);
% xx = abs(xx);
% dx = xx(:).*cos(tt(:));
% dy = xx(:).*sin(tt(:));
% X = ((xx(:) - xx(:)').^2) / sigma_r^2 + ((tt(:) - tt(:)').^2) / sigma_t^2;
% MX = 1/(2*pi*sigma_r*sigma_t)*exp(-X/2) - 1/(2*pi*K^2*sigma_r*sigma_t)*exp(-X/(2*K^2));
% MX = MX./max(MX,[],'all');
% 
% maxIteration = 10;
% %organize population responses
% for indResp = 2:maxIteration
%     %Apply my Weights
%     if indResp <= 2
%         pop_resp{indResp} = reshape(W(:,:,indResp)*reshape(pop_resp{indResp-1},sze(1)*sze(2),[]),sze);
%     %iterates mexican hat weigthing function
%     else
%         CT = reshape(pop_resp{indResp-1},sze(1)*sze(2),[]);
%         CT = CT./max(CT,[],1);
%         th = 1e-2;
%         CT(CT<th) = 0;
%         CT = MX*CT;
%         CT = (CT./max(CT,[],1));
%         %Thresholding
%         CT = 1 ./(1+exp(-slopeLogistic*(CT-centerLogistic)));
%         pop_resp{indResp} = reshape(CT,sze);
%     end
% end
% 
% %% visualize pop activity
% iterNum = 10;
% tmp2 = reshape(pop_resp{iterNum},sze(1),sze(2),sze(3)*sze(4),sze(5),sze(6));
% contrastNum = 3;
% tmp2 = tmp2(:,:,:,:,contrastNum);
% ind_g = 32; %symmetrical oriented gratings
% for j=nVel-2:nVel-2
%     subPopResp3 = tmp2(:,:,ind_g,j);
%     figure, plot_pop_response(subPopResp3,0,0,param.prefVel)
%     colorbar
%     M = sum(sum(subPopResp3));
%     vx = sum(sum(subPopResp3.*(rho(1:8,:).*cos(theta(1:8,:)))))/M;
%     vy = sum(sum(subPopResp3.*(rho(1:8,:).*sin(theta(1:8,:)))))/M;
%     hold on, plot(vx/2,vy/2,'k*')
% end
% %% visualize percetange orientation estimation error with contrast difference increasing
% truetheta = 0;
% % theta_pop = 0:pi/param.nOrient:pi-pi/param.nOrient;
% [theta1,theta2] = meshgrid(linspace(-3*pi/8,3*pi/8,7));
% theta_g = [theta1(:),theta2(:)];
% theta_g(1:8:end,:) = [];
% theta_g1 = reshape(theta_g(:,1),sze(3),sze(4));
% theta_g2 = reshape(theta_g(:,2),sze(3),sze(4));
% %sample points of population response
% x = linspace(min(param.prefVel),max(param.prefVel),101);
% [Xout, Yout] = meshgrid(x);
% thetaout = atan(Yout./Xout);
% thetaout = thetaout + (thetaout<0)*pi;
% thetaout(isnan(thetaout)) = 0;
% rhoout = Xout.*cos(thetaout) + Yout.*sin(thetaout);
% figure,
% for indW=2:maxIteration
% %         for sel = 1:numel(lambda)
%     %select pop_response for one one plaid velocity
%     subPopResp1 = reshape(pop_resp{indW},sze(1),sze(2),sze(3)*sze(4),sze(5),sze(6));
%     subPopResp1 = squeeze(subPopResp1(:,:,:,:,contrastNum));
% %     figure, 
%     ind = 1:49;
%     ind(1:8:49) = [];
%     for i=ind_g:ind_g
%         subPopResp2 = squeeze(subPopResp1(:,:,i,:));
%         for indVel=1:nVel
%             subPopResp3 = subPopResp2(:,:,indVel);
%             M = sum(sum(subPopResp3));
%             vx(indVel) = sum(sum(subPopResp3.*(rho(1:8,:).*cos(theta(1:8,:)))))/M;
%             vy(indVel) = sum(sum(subPopResp3.*(rho(1:8,:).*sin(theta(1:8,:)))))/M;
%             %interpolo per avere una popolazione più fitta
%             tmp = interp2(rho,theta,[subPopResp3; fliplr(subPopResp3(1,:))],rhoout,thetaout);
% % 
% %                 [i,indW]
% %                 pause(0.15)
% %                 close
% % %             %salvami, per ogni valore di velPlaid, l'orientamento con
% % %             %la risp max
% % %             [m,tmpindOr] = max(subPopResp3Interpolated); 
% % %             [m,tmpindVel] = max(m);
% % %             indOr(indVel) = tmpindOr(tmpindVel);
% % %             indVel(indVel) = tmpindVel;
%         end
% %         thetapop = mod(atan2(x(indOr),x(indVel)),2*pi);
%         thetapop = atan(vy./vx) + pi*(param.prefVel(1:nVel)<0);
%         newtruetheta = repmat(truetheta,[1 5]) + pi*(param.prefVel(1:nVel)<0);
%         hold on, plot(param.prefVel(1:nVel), abs(rad2deg( newtruetheta - thetapop )))
%         title(['lambda',num2str(lambda(eta))])
% %         subplot(7,7,ind(i)), plot(param.prefVel(1:nVel), abs(rad2deg( newtruetheta - thetapop )))
% %         ylim([0 50])
% %         xlim([-2 0])
%     end
% end
% legend(strtrim(cellstr(num2str((1:maxIteration-1)'))'))
% end
% rad2deg(truetheta)
% 
% %     pause
% %     close all