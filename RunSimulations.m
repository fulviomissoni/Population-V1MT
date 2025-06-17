%Run simulations 

clear
close all
clc

addpath FUNCTIONS

% NOTE: stim.dur is something that I used to initialise the stimulus
% duration of the used stimulus 
% As regarding the plaid (I used the script initPlaidStimulus) the dur is
% set by the function at 43; for the RDS I will set at the same length;
% !!in any case it should at least match the
% length of the temporal profile of the cells
%%  POPULATION INIT

%SPATIAL FILTERS

% % FILTER 11x11
% filter_file='FILTERS/Gt11B0.0833f0.25.mat';     %Spatial domain - components for 8 orientations of 11x11 Gabor 
% % filter_file='FILTERS/Gt15B0.0833f0.250.mat';     %Spatial domain - components for 8 orientations of 11x11 Gabor 
% k0=0.25;        %SPATIAL FREQUENCY
% n_orient = 8;
% samples = 11;
% filter_sample = 11;
% % samples = 15;
% % RELATIVE BANDWIDTH => B=0.0833;

%FILTER 43x43
filter_file = 'FILTERS/Gt43B0.0208f0.063.mat';
k0 = 0.063;             %SPATIAL FREQUENCY  [cycle/pix]
samples = 127;          %STIMULUS DIMENSION [pix] 
%choice size big enough to obtain good tuning curves in response to RDS
n_orient = 8;
filter_sample = 43;
% RELATIVE BANDWIDTH => B=0.0208;

%TEMPORAL FILTER
ft_choice = 'gabor'; % 'gabor'; 'exp_decay'; 'adelson_bergen'
%PREFERRED VELOCITY
% v = 0;   
v  = linspace(-1,1,11)*2;
% kk = [-3 -1.5  0.5  1.5 3]; %Preferred velocity with Adelson_Bergen

%NORMALIZATION VALUES
% alpha = [1,0.2];

sigma_pool = 3;
num_or_ch_pooled = [8,1];

%Organize the input parameters for the functions
param.spat_freq     = k0;               
param.n_orient      = n_orient;       
param.pref_vel      = v;              
param.temp_filt     = ft_choice;      
param.spatial_filt  = filter_file;    
param.samples       = samples;        
param.filter_sample = filter_sample;
% param.norm_param    = alpha;  
param.sigma_pool    = sigma_pool;
% param.num_or_ch_pooled = num_or_ch_pooled;

param.diff_c = linspace(0,1,11);
diff_c = param.diff_c;
vplaid = 1.8;
truetheta = pi/6 + pi/2;
theta_grat1 = truetheta + pi/4;
theta_grat2 = truetheta + deg2rad(linspace(-45,75,3));

[vgrat(:,:,1), vgrat(:,:,2)] = meshgrid(theta_grat1,theta_grat2);
vgrat = reshape(vgrat,[],2);


dt = datetime('now');
str = char(dt, 'yyyyMMdd_HHmmss');  % Example: "20230603_153045"
%% POP RESPONSE TO TYPE II PLAID 
% This part of the script just simulate the response of normalised V1
% motion-detectors with different contrast levels plaid of type II 
% -) Here we wanted to analyse the effect of normalisation as a function of
% contrast level
% diff_c is used to to simulated the different contrast levels
% a2 to simulate the parameter in the formula of normalisation (Remember
% was: C1 / (a1 + a2/mean(activity) * Sum of activity)
% Sum of activity is modulated with num_or_ch_pooled to modulate the
% normalisation selectivity - Here we just simulate from 8 (all population
% = no selectivity) to 1 (extreme level)
% MT selectivity is obtained with weights which we defined in different
% ways -> see script MT_cells_from_Simulations to see how it works


%normalisation parameters
%Note: I don't work on alpha1 in the normalisation stage (C1/(a1 + a2*pool))

alpha2 = linspace(0,1,11);
alpha1 = [1,zeros(1,10)];
% alpha2 = 1;
% alpha1 = 0;
 
[num_or_pool,a2] = meshgrid(num_or_ch_pooled,alpha2);
% contr = contr(:);
num_or_pool = num_or_pool(:);
a2 = a2(:);
a1 = zeros(length(a2),1);
a1(a2 == 0) = 1;
% num_or_ch_pooled = param.num_or_ch_pooled;
% num_or_pool = num_or_ch_pooled;

%initialise the image
totsim = cell(3,numel(num_or_pool));

for i=1:numel(num_or_pool)
    tic
%      diff_c = contr(i); %contrast difference between gratings

    stim = initPlaidStimulus(truetheta,[vgrat(:,1), vgrat(:,2)],vplaid,diff_c(:),k0,filter_sample,0);
    % stim(:).type = "plaid";
    % stim(:).mode = 1; %implementation mode (see GeneratePlaid)
    % stim.disp = 0;
    % stim.k_gauss = 0;

    %SIMULATION 
    param.num_or_ch_pooled = num_or_pool(i);

    param.norm_param = [a1(i), a2(i)];
    [e,param] = motionPopV1MT(param,stim); % remember 'e' (population activity) will be a matrix of n_complex_cell X n_orient X n_vel X n_stim_parameters (in this case the length of diff_c)
    th = 2e-2;    
    sze_e = size(e);
    
    %DISPLAY RESULTS
    theta_cell_OUT = 0:pi/param.n_orient:pi-pi/param.n_orient;    
    [xx,tt] = meshgrid(param.pref_vel,theta_cell_OUT);
    
    mypath = 'SIMULATIONS\PlaidAnalysis\';
    namesimtmp = [mypath,'tmpSimulationTot_NormEffect_',num2str(i),'_',str];
    save(namesimtmp,'e','stim','param')
    totsim{1,i} = e;
    totsim{2,i} = stim;
    totsim{3,i} = param;
    toc
    disp(['Simulation ',num2str(i), ' finished'])
end
param.norm_param = [alpha1, alpha2];
param.num_or_ch_pooled = num_or_ch_pooled;
mypath = 'SIMULATIONS\PlaidAnalysis\';
namesim = [mypath,'SimulationTot_NormEffect',str];

save(namesim,'totsim')

%% POP RESPONSE to PLAID II with some noise on it
% Here we wanted to simulate what happen in the real-world. In the
% real-word usually stimuli are not so sharp in frequency domain but there
% is some noise (perceptual noise, noise in the background etc..)
% So, we simulate response to plaid II with some noise on it.
% To obtain the cell response we did the average on different locations
% So I used a type II plaid of size 6 times the RF of the neuron and take the
% average 
% % % TO DO: EDIT MOTIONV1MT THE PART OF STIMULUS SELECTION

%normalisation parameters
%Note: I don't work on alpha1 in the normalisation stage (C1/(a1 + a2*pool))

alpha2 = linspace(0,1,11);
alpha1 = [1,zeros(1,10)];
% alpha2 = 1;
% alpha1 = 0;
 
[num_or_pool,a2] = meshgrid(num_or_ch_pooled,alpha2);
% contr = contr(:);
num_or_pool = num_or_pool(:);
a2 = a2(:);
a1 = zeros(length(a2),1);
a1(a2 == 0) = 1;
% num_or_ch_pooled = param.num_or_ch_pooled;
% num_or_pool = num_or_ch_pooled;

%initialise the image
totsim = cell(3,numel(num_or_pool));

for i=1:numel(num_or_pool)
    tic
%      diff_c = contr(i); %contrast difference between gratings

    stim = initPlaidStimulus(truetheta,[vgrat(:,1), vgrat(:,2)], ...
        vplaid,diff_c(:),k0,filter_sample,0 ...
        ,'type',"plaid_noise", ...
        'sigma_noise',1);
    % stim(:).type = "plaid";
    % stim(:).mode = 1; %implementation mode (see GeneratePlaid)
    % stim.disp = 0;
    % stim.k_gauss = 0;

    %SIMULATION 
    param.num_or_ch_pooled = num_or_pool(i);

    param.norm_param = [a1(i), a2(i)];
    [e,param] = motionPopV1MT(param,stim); % remember 'e' (population activity) will be a matrix of n_complex_cell X n_orient X n_vel X n_stim_parameters (in this case the length of diff_c)
    th = 2e-2;    
    sze_e = size(e);
    
    %DISPLAY RESULTS
    theta_cell_OUT = 0:pi/param.n_orient:pi-pi/param.n_orient;    
    [xx,tt] = meshgrid(param.pref_vel,theta_cell_OUT);
    
    mypath = 'SIMULATIONS\PlaidAnalysis\withNoise\';
    namesimtmp = [mypath,'tmpSimulationTot_Noise_NormEffect_',num2str(i),'_',str];
    save(namesimtmp,'e','stim','param')
    totsim{1,i} = e;
    totsim{2,i} = stim;
    totsim{3,i} = param;
    toc
    disp(['Simulation ',num2str(i), ' finished'])
end
param.norm_param = [alpha1, alpha2];
param.num_or_ch_pooled = num_or_ch_pooled;
mypath = 'SIMULATIONS\PlaidAnalysis\withNoise\';
namesim = [mypath,'SimulationTot_Noise_NormEffect',str];

save(namesim,'totsim')

%% POP RESPONSE TO MOVING RDS
% This is related to the same of step of before (having a broader frequency
% response) and see what happen to the cells
% % % % TO DO: I don't remember how I did


% Parameters of the simulation
clear stim
stim.type = 'RDS_moving';
stim.disp = 0;
stim.truetheta = 0;
stim.vgrat =  [vplaid , 0]; %velocity in px/sec; 
stim.dur = 43;
stim.theta_g = [0, 0];
num_or_ch_pooled = 1;
% param.diff_c = linspace(0,1,11);
% diff_c = param.diff_c;
%norm parameters
% alpha2 = linspace(1,0,11);
alpha2 = 1;
alpha1 = 0;
[num_or_pool,a2] = meshgrid(num_or_ch_pooled,alpha2);
% contr = contr(:);
num_or_pool = num_or_pool(:);
a2 = a2(:);
% num_or_ch_pooled = param.num_or_ch_pooled;
% num_or_pool = num_or_ch_pooled;

%initialise the image
totsim = cell(3,numel(num_or_pool));

for i=1:numel(num_or_pool)
    tic
    
    %SIMULATION 
    param.num_or_ch_pooled = num_or_pool(i);
    param.norm_param = [alpha1, a2(i)];
    [e,param] = motionPopV1MT(param,stim); % remember 'e' (population activity) will be a matrix of n_complex_cell X n_orient X n_vel X n_stim_parameters (in this case the length of diff_c)
    th = 2e-2;    
    sze_e = size(e);
    
    %DISPLAY RESULTS
    theta_cell_OUT = 0:pi/param.n_orient:pi-pi/param.n_orient;    
    [xx,tt] = meshgrid(param.pref_vel,theta_cell_OUT);
    
    mypath = 'SIMULATIONS\RDSAnalysis\';
    namesimtmp = [mypath,'tmpSimulationTot_NormEffect_',num2str(i),'_',str];
    save(namesimtmp,'e','stim','param')
    totsim{1,i} = e;
    totsim{2,i} = stim;
    totsim{3,i} = param;
    toc
    disp(['Simulation ',num2str(i), ' finished'])
end
param.norm_param = [ repmat(alpha1,1,length(alpha2)), alpha2];
param.num_or_ch_pooled = num_or_ch_pooled;
mypath = 'SIMULATIONS\RDSAnalysis\';
namesim = [mypath,'SimulationTot_NormEffect',str];

save(namesim,'totsim')