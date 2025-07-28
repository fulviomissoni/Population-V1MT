clear variables
close all

addpath FUNCTIONS\
% load("SIMULATIONS\PlaidAnalysis\NoNoise\SimulationTot_NormEffect20250709_122400.mat")
load("SIMULATIONS\PlaidAnalysis\Noise\SimulationTot_Noise_NormEffect20250722_122439.mat")

% M = 3e5; %compute from population response to plaid
% M = max(data,[],"all");
data = cat(5,totsim{1,:});
stim = squeeze(cat(4,totsim{2,:}));
param = squeeze(cat(4,totsim{3,:}));

% num_or_ch_pooledTot = unique(cat(1,param.num_or_ch_pooled));
num_diff_c = unique(cat(1,param.diff_c));
norm_param_sigma = cat(1,param.sigma_orient);
% num_norm = unique(norm_param(:,2));

sze = size(data);
%pop parameters
numComplexCellsPopulations = sze(1);
numOrient = sze(2);
numVel = sze(3);
numTestedPopParameters = sze(5); %num population parameters tested

%stim tested parameters
numTestedStim = sze(4);       %stim parameters tested
TotTestedGrat = cat(1,stim.theta_g); %Theta in radiant
TestedGrat_1 = unique(TotTestedGrat(:,1));
TestedGrat_2 = unique(TotTestedGrat(:,2));
TotTestedVGrat = cat(1,stim.vgrat); %in px/sec
TestedVGrat_1 = unique(TotTestedVGrat(:,1));
TestedVGrat_2 = unique(TotTestedVGrat(:,2));
TotTestedDiff = diff(cat(1,stim.contrast),[],2); 
TestedDiff = unique(TotTestedDiff);

TotTestedGrat = reshape(TotTestedGrat, ...
        [numel(TestedGrat_1),numel(TestedGrat_2) ...
        numel(TestedDiff),numTestedPopParameters,2]);
TotTestedDiff = reshape(TotTestedDiff, ...
        [numel(TestedGrat_1),numel(TestedGrat_2) ...
        numel(TestedDiff),numTestedPopParameters]);

% testedTheta_grat = rad2deg( ...
%     reshape( ...
%     cat(1,stim.theta_g), ...
%     [size(stim,1),numel(num_or_ch_pooledTot),numel(num_diff_c),2]));
% testedV_grat = abs( ...
%     reshape( ...
%     cat(1,stim.vgrat), ...
%     [size(stim,1),numel(num_or_ch_pooledTot),numel(num_diff_c),2]));
% testedContrasts = reshape( ...
%     cat(1,stim.contrast), ...
%     [size(stim,1),numel(num_or_ch_pooledTot),numel(num_diff_c),2]);
% testedTrueThetaPlaid = rad2deg( ...
%     reshape( ...
%     cat(1,stim.truetheta), ...
%     [size(stim,1),numel(num_or_ch_pooledTot),numel(num_diff_c)]));

data_reshaped = reshape(data,[sze(1:end-2), ...
    numel(TestedGrat_1),numel(TestedGrat_2), ...
    numel(TestedDiff),numTestedPopParameters]);
% Just to remember the size of the matrix is this:
% filter_polarity_inPop x num_orient x num_vel x num_tested_grat1 x num_tested_grat2 
% x num_tested_contrast_diff x num_tested_pop_parameters =  size(data_reshaped)
% data = logistic(cat(5,totsim{1,:}),2.5e2,.003);


%population tested parameters
plotTestedVelocitiesDirection

%% Look at the data

% figure,
% stim(1).disp = 1;
% pl = plaid(stim(1));
% generate_plaid(pl);
% Generate the population response (figure) for each combination of parameters
%each figure is the pop_response without normalisation and same contrast
%level of gratings for 3 different combinations of gratings directions
%(column) and 4 different weighting system

% for ii = 1:numel(num_or_ch_pooledTot)
%     for jj = 1:numel(num_norm)
sze_reshaped = size(data_reshaped);
data_reshaped = reshape(data_reshaped,[sze_reshaped(1:3),prod(sze_reshaped(4:5)),sze_reshaped(6:end)]);
theta_cell_OUT = 0:pi/param(1).n_orient:pi-pi/param(1).n_orient;    
[xx,tt] = meshgrid(param(1).pref_vel,theta_cell_OUT);
%Get the different weighting models and repeated for each tested values
W = permute( ...
    repmat( ...
    weightingfunctions(xx,tt), ...
    1,1,1,prod(sze_reshaped(4:5)),sze_reshaped(6),sze_reshaped(7)), ...
    [1,2,4,5,6,3]);  

%for each orientation combination make a figure
%I plot what happen if I have the grating contrast level at different ratio
%(columns) with different normalisation levels (rows)
%EVEN POPULATION
for ii = 1:prod(sze_reshaped(4:5))
    
    e_1(:,:,ii,:,:) = squeeze(data_reshaped(1,:,:,ii,:,:));
    e_2(:,:,ii,:,:) = squeeze(data_reshaped(2,:,:,ii,:,:));

    w = squeeze(W(:,:,ii,:,:,:));
    pop_resp_odd(:,:,ii,:,:,:) = getpopresponse(w,e_1(:,:,ii,:,:)); % Apply weighting to the response data
    pop_resp_even(:,:,ii,:,:,:) = getpopresponse(w,e_2(:,:,ii,:,:)); % Apply weighting to the response data
    
    % pop_resp = squeeze(W .* e); 
    pop_resp = pop_resp_even - pop_resp_odd;
    % figure, popresponse_tiled(e_2(:,:,ii,:,:,:));
    % figure, popresponse_tiled(squeeze(pop_resp_even(:,:,ii,:,:,2)));

    % figure, popresponse_tiled(cat(4,e(:,:,:,1),pop_resp))
end
figure, plot(TestedDiff), hold on,plot(norm_param_sigma)

%% Decode Activity
%Decoding activity of pop_resp with centre of mass of activity
%Before, activity is processed iteratively with mexican hat (defined by
%sigma_r,sigma_t,K parameters) in max_iteration iteration numbers
%each iteration the activity is thresholded with logistic function defined
%by (logistic_slope and logistic_centre parameters)
%
% n_orient x n_vel x = size(pop_resp) 

%Decoding Parameters
%Mexican Hat parameters
%Remember that a mexican hat is defined as 
% MX = 1/(2*pi*sigma_r*sigma_t)*exp(-X/2) - ...
%    1/(2*pi*K^2*sigma_r*sigma_t)*exp(-X/(2*K^2));
% sigma_r = 0.1;  
% sigma_t = 0.4;  
% K = 1.5; %inihibitory field size factor
% % Thresholding parameter
% 
% logistic_slope = 8;
% logistic_centre = 0.65;
% max_iteration = 105;

sigma_r = 0.2;  
sigma_t = 0.3;  
K = 2; %inihibitory field size factor
% Thresholding parameter
 
logistic_slope = 6;
logistic_centre = 0.3;
max_iteration = 12;


sze_pop = size(pop_resp_even);
num_cond = prod(sze_pop(3:end));
%Activity Decoding
for ii = 1:num_cond
    [PR_decoded(:,:,:,ii),vx(ii),vy(ii)] = DecodeMxHat( ...
                                squeeze(pop_resp_even(:,:,ii)), ...  %Activity
                                param(1), ...                                  %Population Parameters
                                sigma_r, ...                                %see code
                                sigma_t, ...
                                K, ...
                                max_iteration, ...
                                logistic_slope, ...
                                logistic_centre);
end

PR_decoded = reshape(PR_decoded,[sze_pop(1),sze_pop(2),max_iteration,sze_pop(3:end)]);
vx = reshape(vx,sze_pop(3:end));
vy = reshape(vy,sze_pop(3:end));

%% plot
figure, popresponse_tiled(squeeze(PR_decoded(:,:,1,1,:,1,1)))
% tiledlayout('flow')
figure,
% for ii = 1:6
% v = squeeze(sqrt(vx(:,:,1,1).^2 + vy(:,:,1,1).^2));
v = squeeze(sqrt(vx.^2 + vy.^2));

plot(TestedDiff,v), legend(string(rad2deg(theta2)))
%legend(string(TestedDiff(1))),
ylim([1,1.8])
% end