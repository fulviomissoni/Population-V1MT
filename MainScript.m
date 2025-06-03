clear
close all
clc

addpath FUNCTIONS
%%  POPULATION INIT

%SPATIAL FILTERS

% % FILTER 11x11
% filter_file='FILTERS/Gt11B0.0833f0.25.mat';     %Spatial domain - components for 8 orientations of 11x11 Gabor 
% % filter_file='FILTERS/Gt15B0.0833f0.250.mat';     %Spatial domain - components for 8 orientations of 11x11 Gabor 
% k0=0.25;        %SPATIAL FREQUENCY
% n_orient = 8;
% samples = 11;
% % samples = 15;
% % RELATIVE BANDWIDTH => B=0.0833;

%FILTER 43x43
filter_file = 'FILTERS/Gt43B0.0208f0.063.mat';
k0 = 0.063;             %SPATIAL FREQUENCY  [cycle/pix]
samples = 127;          %STIMULUS DIMENSION [pix] 
%choice size big enough to obtain good tuning curves in response to RDS
n_orient = 8;
% RELATIVE BANDWIDTH => B=0.0208;

%TEMPORAL FILTER
ft_choice = 'gabor'; % 'gabor'; 'exp_decay'; 'adelson_bergen'
%PREFERRED VELOCITY
% v = 0;   
v  = linspace(-1,1,11)*2;
% kk = [-3 -1.5  0.5  1.5 3]; %Preferred velocity with Adelson_Bergen

%NORMALIZATION VALUES
alpha = [1,0.2];
sigma_pool = 3;
num_or_ch_pooled = [8,1];

%Organize the input parameters for the functions
param.spat_freq     = k0;               
param.n_orient      = n_orient;       
param.pref_vel      = v;              
param.temp_filt     = ft_choice;      
param.spatial_filt  = filter_file;    
param.samples       = samples;        
param.norm_param    = alpha;  
param.sigma_pool    = sigma_pool;
param.num_or_ch_pooled = num_or_ch_pooled;

%% POP RESPONSE TO TYPE II PLAID
param.diff_c = 0:0.1:1;
diff_c = param.diff_c;
[contr,num_or_pool] = meshgrid(diff_c,num_or_ch_pooled);
contr = contr(:);
num_or_pool = num_or_pool(:);
num_or_ch_pooled = param.num_or_ch_pooled;
% num_or_pool = num_or_ch_pooled;
for i=1:numel(num_or_pool)
%      diff_c = contr(i); %contrast difference between gratings
    stim = initStimulus(pi/4,[-3*pi/8, -pi/4],1.8,contr(i));
    stim.mode = 1;
    stim.disp = 0;
    
    %SIMULATION 
    param.num_or_ch_pooled = num_or_pool(i);
    [e,param] = motionPopV1MT(param,stim);
    th = 2e-2;    
    sze_e = size(e);
    
    %DISPLAY RESULTS
    theta_cell_OUT = 0:pi/param.n_orient:pi-pi/param.n_orient;    
    [xx,tt] = meshgrid(param.pref_vel,theta_cell_OUT);
    
    %load tuning curves
%     load('SIMULATIONS/vel_tuning_polarRDS.mat','W')
%     W = W - W.*eye(size(W));
% %     W = W/max(W(:));
%     sigmaGabor = 0.5*4/3*xx(:)'; %sigma = 4/3*xx(:)';
%     sigma = 0.25;
    % W2 = exp(-(xx(:).*cos(tt(:)-tt(:)') - xx(:)').^2./(2*sigma.^2)); %Andre pesi originali (inviluppo exp(.))
    % W2 = exp(-(xx).^2/(2*sigma^2)).*cos(tt(:)-tt(:)'); %Chessa-Solari model
    % W2 = cos(xx(:).*cos(tt(:)'-tt(:)) - xx(:)'); %primo modello per introdurre l'opponenza (inviluppo coseno)
    W2 = exp(-(xx(:).*cos(tt(:)'-tt(:)) - xx(:)').^2./(2*sigmaGabor.^2)).* ...
    cos(2*pi*1./(4*xx(:)').*(xx(:).*cos(tt(:)'-tt(:)) - xx(:)')); %inviluppo Gabor per gestire la quantità e la posizione dei pesi negativi
   % W2 = reshape(reshape(W2,8,11,8,11)./max(reshape(W2,8,11,8,11),[],4),88,88);
    W2(isnan(W2)) = 0;
    
    pop_resp(:,:,i) = squeeze(e(2,:,:));
    sze = size(squeeze(pop_resp(:,:,i)));
    pop_resp_TunCurves(:,:,i) = reshape((W*reshape(squeeze(pop_resp(:,:,i)),sze(1)*sze(2),[])),sze);
    pop_resp_BioGautama(:,:,i) = reshape((W2*reshape(squeeze(pop_resp(:,:,i)),sze(1)*sze(2),[])),sze);
end
param.num_or_ch_pooled = num_or_ch_pooled;
OldFolder = cd;
cd('SIMULATIONS\PlaidAnalysis')
save('HIGHSIZEpop_resp_a2_0_2','pop_resp','pop_resp_BioGautama','pop_resp_TunCurves','W2','stim','param')
cd(OldFolder)
%% PLAID VELOCITY - POPULATION ENCODING

% Gaussian extension parameters
sigma_r = 0.1;  
% sigma_t = 0.4;  
sigma_t = 0.4;
K = 1.5; %inihibitory field size factor
% Thresholding parameter
logistic_slope = 9; 
logistic_centre = 0.7;
max_iteration = 6;

[pop_resp_V1MT,vx,vy] = DecodeMxHat(pop_resp_BioGautama,param,sigma_r,sigma_t,K,max_iteration,logistic_slope,logistic_centre);
% MEXICAN HAT WEIGHTS
truetheta = repmat(stim.truetheta,[numel(diff_c), max_iteration, numel(num_or_ch_pooled)]);
figure
for i=1:max_iteration
    subplot(1,max_iteration,i)
    plot(diff_c,abs(rad2deg(squeeze(truetheta(:,i,:)-atan2(vy(:,i,:),vx(:,i,:))))))
    ylim([0,max(max(max(abs(rad2deg(squeeze(truetheta-atan2(vy,vx)))))))])
    if i==1
        ylabel('\Delta\theta')
    end
    if i==ceil(max_iteration/2)
        xlabel('\Deltac')
    end
    if i==max_iteration
        legend(['numOrPooled = ',num2str(param.num_or_ch_pooled(1))],['numOrPooled = ',num2str(param.num_or_ch_pooled(2))])
    end
%     pause
end
%SAVE
% saveas(gcf,'angle error plot','fig')
% 
figure, h=plot(diff_c,abs(rad2deg(squeeze(truetheta(:,:,1)-atan2(vy(:,:,1),vx(:,:,1))))));
title('numOrPooled = 1')
figure, plot(diff_c,abs(rad2deg(squeeze(truetheta(:,:,2)-atan2(vy(:,:,2),vx(:,:,2))))))
title('numOrPooled = 8')

%% EXAMPLE: POP ACTIVITY ON SET OF TYPE II PLAID (WITH DIFFERENTS VEL VECTOR; ALL COMBINATIONS)

%STIMULUS DEFINITION
truetheta = 0;
plaid_vel = param.pref_vel(1:5);
[theta1, theta2] = meshgrid(linspace(-3*pi/8,3*pi/8,7));
theta_g = [theta1(:), theta2(:)];
theta_g(1:8:end,:) = [];
%differences in contrast sensitivity
diff_contrast = 0:0.2:0.6;
stim = initStimulus(truetheta(:),theta_g,plaid_vel,diff_contrast);

stim.dur = 47; %duration in frame
stim.mode = 1;
stim.disp = 0; %set to 1 to show visual stimulus in a figure
%% SIMULATION

lambda = [0,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,1e1,1e2];
for j = 1:numel(lambda)
    param.norm_param = [1;lambda(j)];
    [e, param] = motionPopV1MT(param,stim);
    th = 2e-2;
    %SAVE DATA
    path = 'SIMULATIONS';
    OldFolder = cd;
    cd(path);
    filename = ['Vel_tuning_PlaidII_lambda',num2str(lambda(j)),'_diffContrasts'];
    save(filename,'e','param','stim','-v7.3')
    cd(OldFolder)
end

figure, plot(diff_c,abs(rad2deg(repmat(stim.truetheta,[numel(diff_c) numel(num_or_ch_pooled)])-atan2(vy,vx))))
%%  LOCAL FUNCTIONS
%initialize stimulus parameters

%% reshape stimulus function
%reshape parameters, from ngrid to n-d matrices