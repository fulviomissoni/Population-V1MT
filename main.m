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
% samples = 11;
% % samples = 15;
% % RELATIVE BANDWIDTH => B=0.0833;

%FILTER 43x43
filter_file = 'FILTERS/Gt43B0.0208f0.063.mat';
k0 = 0.063;             %SPATIAL FREQUENCY  [cycle/pix]
samples = 420;          %STIMULUS DIMENSION [pix]
% RELATIVE BANDWIDTH => B=0.0208;

%PHASE SHIFT
% ph_ = pi/4:-pi/4:-5/4*pi; %this choice is related to the offset value of phase_shift (pi/2)
ph_ = pi/2;
n_orient = 8;
ph = repmat(ph_,n_orient,1);  %phase_shift value for each orientation channel 

%OCULAR DOMINANCE
a = 1; 
%TEMPORAL FILTER
Ft_choice = 'gabor'; % 'gabor'; 'exp_decay'; 'adelson_bergen'
%PREFERRED VELOCITY
% v = 0;   
v  = linspace(-1,1,11)*2;
% kk = [-3 -1.5  0.5  1.5 3]; %Preferred velocity with Adelson_Bergen

%NORMALIZATION VALUES
alpha = [1;1e-4];

%Organize the input parameters for the functions
param.spatFreq    = k0;             
param.ocDom       = a;              
param.phShift     = ph;             
param.nOrient     = n_orient;       
param.prefVel     = v;              
param.tempFilt    = Ft_choice;      
param.spatialFilt = filter_file;    
param.samples     = samples;        
param.normParam   = alpha;          

%% POP RESPONSE TO TYPE II PLAID
stim = init_stimulus(0,[3*pi/8 -3*pi/8],param.prefVel(end),0.4);
stim.mode = 1;
stim.disp = 0;

%SIMULATION
[e,param] = motion_popV1MT(param,stim);
th = 2e-2;

%DISPLAY RESULTS
theta_cell_OUT = 0:pi/param.nOrient:pi-pi/param.nOrient;

[xx,tt] = meshgrid(param.prefVel,theta_cell_OUT);

% load 'SIMULATIONS/GautamaWeights88_Plaid.mat'
%Explicit intersection of constraints method to compute weigths
W2 = exp(-(xx(:).*cos(tt(:)'-tt(:)) - xx(:)').^2/(2*0.25^2));
W2 = W2 - eye(size(W2));

pop_resp = squeeze(e(3,66/2,66/2,:,:));
sze = size(pop_resp);
pop_resp_BioGautama = reshape(W*reshape(pop_resp,sze(1)*sze(2),[]),sze);
pop_resp_BioGautama2 = reshape((W2*reshape(pop_resp,sze(1)*sze(2),[])),sze);

%POP RESPONSE TO MOVING RDS (HORIZONTALLY AT 2 pix/frame)
figure,plot_pop_response(pop_resp,0,0,param.prefVel)
title('POP RESPONSE')
figure,plot_pop_response(pop_resp_BioGautama,0,0,param.prefVel)
title('POP RESP BIO GAUTAMA')
figure,plot_pop_response(pop_resp_BioGautama2,0,0,param.prefVel)
title('POP RESP BIO GAUTAMA2')

%% EXAMPLE: POP ACTIVITY ON SET OF TYPE II PLAID (WITH DIFFERENTS VEL VECTOR; ALL COMBINATIONS)

%STIMULUS DEFINITION
truetheta = 0;
plaid_vel = param.prefVel(1:5);
[theta1,theta2] = meshgrid(linspace(-3*pi/8,3*pi/8,7));
theta_g = [theta1(:),theta2(:)];
theta_g(1:8:end,:) = [];
%differences in contrast sensitivity
diff_contrast = 0:0.2:0.6;
stim = init_stimulus(truetheta(:),theta_g,plaid_vel,diff_contrast);

stim.dur = 47; %duration in frame
stim.mode = 1;
stim.disp = 0; %set to 1 to show visual stimulus in a figure
%% SIMULATION

lambda = [0,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,1e1,1e2];
for i = 1:numel(lambda)
    param.normParam = [1;lambda(i)];
    [e,param] = motion_popV1MT(param,stim);
    th = 2e-2;
    %SAVE DATA
    path = 'SIMULATIONS';
    OldFolder = cd;
    cd(path);
    filename = ['newVel_tuning_PlaidII_lambda',num2str(lambda(i)),'_difContrasts'];
    save(filename,'e','param','stim','-v7.3')
    cd(OldFolder)
end

%%  FUNCTIONS
function stim = init_stimulus(truetheta,theta,vpld,diff_contrast)
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