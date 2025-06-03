function [pop_resp_V1MT, vx, vy] = DecodeMxHat(pop_resp,param,sigma_r,sigma_t,K,max_iteration,logistic_slope,logistic_centre)
%
%Decoding activity of pop_resp_BioGautama with centre of mass of activity
%Before, activity is processed iteratively with mexican hat (defined by
%sigma_r,sigma_t,K parameters) in max_iteration iteration numbers
%each iteration the activity is thresholded with logistic function defined
%by (logistic_slope and logistic_centre parameters)
%
% 
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
