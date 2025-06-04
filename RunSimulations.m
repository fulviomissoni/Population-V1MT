%Run simulations 

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
dt = datetime('now');
str = char(dt, 'yyyyMMdd_HHmmss');  % Example: "20230603_153045"
namesim = ['HIGHSIZEpop_resp_a2_0_2',str];
save(namesim,'pop_resp','pop_resp_BioGautama','pop_resp_TunCurves','W2','stim','param')
cd(OldFolder)