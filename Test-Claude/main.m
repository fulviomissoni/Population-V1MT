
%Run simulations 

clear
close all
clc

% addpath FUNCTIONS
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
a2 = linspace(0,1,10);

%Organize the input parameters for the functions
param.spat_freq     = k0;               
param.n_orient      = n_orient;       
param.pref_vel      = v;              
param.temp_filt     = ft_choice;      
param.spatial_filt  = filter_file;    
param.samples       = samples;        
param.norm_param    = [ones(1,10); a2];  
param.sigma_pool    = sigma_pool;
param.num_or_ch_pooled = n_orient;

%% POP RESPONSE TO TYPE II PLAID
param.diff_c = 0:0.1:1;
diff_c = param.diff_c;
[contr,alpha2] = meshgrid(diff_c,a2);
contr = contr(:);
alpha2 = alpha2(:);
vpld = 1.8; %deg/sec (?)

%STIM
stim = initStimulus(pi/4,[-3*pi/8, -pi/4],vpld,contr);
stim.mode = 1;
stim.disp = 0;

II = makePlaids()

results = run_parameter_sweep_vectorized(II, your_base_param, ...
    'contr', contr, 'norm_alpha', alpha2);

% Analisi automatica risultati
analyze_results(results, 'contr', param1_values, 'norm_alpha', param2_values);
