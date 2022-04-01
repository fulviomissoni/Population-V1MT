clear
% close all
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
samples = 420;          %STIMULUS DIMENSION [pix] 
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
alpha = [1;1];
sigma_pool = 3;
num_or_ch_pooled = 1;

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
stim = initStimulus(0,[3*pi/8 -3*pi/8],param.pref_vel(end),0.2);
stim.mode = 1;
stim.disp = 0;

%SIMULATION
[e,param] = motionPopV1MT(param,stim);
th = 2e-2;

sze_e = size(e);
%DISPLAY RESULTS
theta_cell_OUT = 0:pi/param.n_orient:pi-pi/param.n_orient;

[xx,tt] = meshgrid(param.pref_vel,theta_cell_OUT);

%Explicit intersection of constraints method to compute weigths
W2 = exp(-(xx(:).*cos(tt(:)'-tt(:)) - xx(:)').^2/(2*0.25^2));
W2 = W2 - eye(size(W2));

pop_resp = squeeze(e(2,:,:));
sze = size(pop_resp);
pop_resp_BioGautama = reshape((W2*reshape(pop_resp,sze(1)*sze(2),[])),sze);

%POP RESPONSE TO MOVING RDS (HORIZONTALLY AT 2 pix/frame)
figure, plotPopResponse(pop_resp,0,0,param.pref_vel)
title('POP RESPONSE')
figure, plotPopResponse(pop_resp_BioGautama,0,0,param.pref_vel)
title('POP RESP BIO GAUTAMA')

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
for i = 1:numel(lambda)
    param.norm_param = [1;lambda(i)];
    [e, param] = motionPopV1MT(param,stim);
    th = 2e-2;
    %SAVE DATA
    path = 'SIMULATIONS';
    OldFolder = cd;
    cd(path);
    filename = ['Vel_tuning_PlaidII_lambda',num2str(lambda(i)),'_diffContrasts'];
    save(filename,'e','param','stim','-v7.3')
    cd(OldFolder)
end

%% PLAID VELOCITY - POPULATION ENCODING
% MEXICAN HAT WEIGHTS
% Gaussian extension parameters
sigma_r = 0.1;  
sigma_t = 0.4;  
K = 1.5; %inihibitory field size factor
% Thresholding parameter
logistic_slope = 9;  
logistic_center = 0.7;

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
max_iteration = 10;

tmp = pop_resp_BioGautama;
%organize population responses
pop_resp_V1MT = pop_resp_BioGautama;
for indResp = 2:max_iteration
    %Apply my Weights
    %iterates mexican hat weigthing function
    CT = reshape(tmp,sze(1)*sze(2),[]);
    CT = CT./max(CT,[],1);
    CT(CT<th) = 0;
    CT = MX*CT;
    CT = (CT./max(CT,[],1));
    %Thresholding
    CT = 1 ./(1+exp(-logistic_slope*(CT-logistic_center)));
    pop_resp_V1MT(:,:,indResp) = reshape(CT,sze);
    tmp = CT;
end

%% VISUALIZE POP ACITIVITY
[rho, theta] = meshgrid(param.pref_vel,linspace(0,pi,n_orient+1));

for j=1:max_iteration
    subPopResp = pop_resp_V1MT(:,:,j);
    figure, plotPopResponse(subPopResp,0,0,param.pref_vel)
    colorbar
    
    M = sum(sum(subPopResp));
    %centre of mass decoding
    vx = sum(sum(subPopResp.*(rho(1:8,:).*cos(theta(1:8,:)))))/M;
    vy = sum(sum(subPopResp.*(rho(1:8,:).*sin(theta(1:8,:)))))/M;
    hold on, plot(vx/2,vy/2,'k*')
end

%%  LOCAL FUNCTIONS
function stim = initStimulus(truetheta,theta,vpld,diff_contrast)
    stim.type = 'plaid';
    stim.truetheta =  truetheta(:); %true orientation
    stim.theta_g = [theta]; %true orientation
    stim.vel_stim = [vpld(:)];
    [x, y, z, c1] = ndgrid(stim.truetheta, stim.theta_g(:,1), stim.vel_stim, 0.5-diff_contrast/2); 
    y = pagetranspose(y);   %pagetranspose serve per rendere le grid create con ndgrid nello stesso ordine dei meshgrid
    c1 = pagetranspose(c1);
    [stim.truetheta, y2, stim.vel_stim,c2] = ndgrid(stim.truetheta, stim.theta_g(:,2), stim.vel_stim, 0.5+diff_contrast/2);
    stim.truetheta = pagetranspose(stim.truetheta);
    y2 = pagetranspose(y2);
    stim.vel_stim = pagetranspose(stim.vel_stim);
    c2 = pagetranspose(c2);
    stim.truetheta = stim.truetheta(:); 
    stim.theta_g = stim.truetheta + [y(:) y2(:)];
    stim.contrast_g = [c1(:),c2(:)];
    stim.vel_stim = stim.vel_stim(:);
    if size(stim.theta_g,2)==2
        stim.ori = [cos(stim.theta_g(:,1)-stim.truetheta), cos(stim.theta_g(:,2)-stim.truetheta)];
        stim.vgrat = [stim.ori(:,1).*stim.vel_stim, stim.ori(:,2).*stim.vel_stim];
        stim.vgrat = round(stim.vgrat,5,"decimals");
    else
        stim.vgrat = stim.vel_stim;
    end
    stim.dur = 43; %duration in frame
    stim.mode = 1;
    stim.disp = 0; %set to 1 to show visual stimulus in a figure
end