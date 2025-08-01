clear variables
close all

addpath FUNCTIONS\
% load("SIMULATIONS\PlaidAnalysis\NoNoise\SimulationTot_NormEffect20250722_122439.mat")
% load("SIMULATIONS\PlaidAnalysis\NoNoise\SimulationTot_NormEffect20250730_183038.mat")
% load("SIMULATIONS\PlaidAnalysis\NoNoise\SimulationTot_NormEffect20250731_142016.mat")
% load("SIMULATIONS\PlaidAnalysis\NoNoise\SimulationTot_NormEffect20250731_164043.mat") %TO DELETE
% load("SIMULATIONS\PlaidAnalysis\NoNoise\SimulationTot_NormEffect20250731_164527.mat")
load("SIMULATIONS\PlaidAnalysis\Noise\SimulationTot_Noise_NormEffect20250801_132719.mat")

% load("SIMULATIONS\PlaidAnalysis\Noise\SimulationTot_Noise_NormEffect20250722_122439.mat")
% load("SIMULATIONS\PlaidAnalysis\Noise\SimulationTot_Noise_NormEffect20250730_175259.mat")
% load("SIMULATIONS\PlaidAnalysis\Noise\SimulationTot_Noise_NormEffect20250731_142016.mat")
% load("SIMULATIONS\PlaidAnalysis\Noise\SimulationTot_Noise_NormEffect20250731_164527.mat") %TO DELETE

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
TotTestedVGrat = reshape(TotTestedVGrat, ...
        [numel(TestedVGrat_1),numel(TestedVGrat_2) ...
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

testedTrueThetaPlaid = unique(stim(1).truetheta);
testedVplaid = unique(stim(1).vpld);
%population tested parameters
% plotTestedVelocitiesDirection
plotPlaidVectors(repmat(testedTrueThetaPlaid,numel(TestedGrat_2),1),...
            squeeze(TotTestedGrat(1,:,1,1,:)),...
            repmat(testedVplaid,numel(TestedGrat_2),1), ...
            squeeze(TotTestedVGrat(1,:,1,1,:))) 
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
    % figure, popresponse_tiled(e_2);
    % figure, popresponse_tiled(squeeze(pop_resp_even(:,:,:,:,2)));

    % figure, popresponse_tiled(cat(4,e(:,:,:,1),pop_resp))
end
% figure, plot(TestedDiff), hold on,plot(norm_param_sigma)

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

% pop_resp = pop_resp_even;
% If your data grid is, say, 8x11:
sigma_r = (param(1).pref_vel(end) - param(1).pref_vel(end-1));
sigma_t = 5/4*(pi/8);
K = 1.8;
logistic_slope = 5;
logistic_centre = 0.6;
max_iteration = 5;


sze_pop = size(pop_resp_even);
num_cond = prod(sze_pop(3:end));
weight_mode = 2;
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

figure, t = popresponse_tiled(cat(3,reshape(squeeze(e_2(:,:,7,:,1)),[sze_pop(1:2),1,sze_pop(4)]),...
                    squeeze(PR_decoded(:,:,:,7,:,1,weight_mode)))); 
title(t,'Population Activity Decoded (rows are different Iteration, col are contrast)')
vx = reshape(vx,sze_pop(3:end));
vy = reshape(vy,sze_pop(3:end));
gradientMap = jet(numel(norm_param_sigma));

n_orient = 8;    
x = linspace(1,n_orient,100);
x_8 = round(linspace(1,100,n_orient));

figure,
w_o = exp(-(x(x_8(3))-(1:n_orient)).^2./(2.*norm_param_sigma.^2));
plot(w_o'),colororder(gradientMap), legend(string(norm_param_sigma)),
title('Tested Orientation Pooling (sample for one neuron)')

v = sqrt(vx.^2 + vy.^2);
ori = atan2(vy,vx);
% figure, plot(squeeze(v(2,:,:,2)))
for ii = 1:size(ori,1)
    figure
    subplot(2,1,1),
    hold on,
    estim = squeeze(ori(ii,:,:,weight_mode)); 
    true = repmat(stim(1,1).truetheta,size(estim,1),size(estim,2));    
    scatter(TestedDiff, abs(rad2deg(angdiff(estim,true))),'filled', ...
        'MarkerFaceAlpha',0.4,'MarkerEdgeColor','none','Color',gradientMap), %legend(string(rad2deg(TestedGrat_2)))
    plot(TestedDiff, abs(rad2deg(angdiff(estim,true)))); legend(string(norm_param_sigma))
    colororder(gradientMap)
    ylim([0,150]),grid on
    title('angle error')

    %velocity
    estim = squeeze(v(ii,:,:,weight_mode)); 
    true = repmat(stim(1,1).vpld,size(estim,1),size(estim,2));
    subplot(2,1,2),
    hold on,
    scatter(TestedDiff, abs(true-estim),'filled', ...
        'MarkerFaceAlpha',0.4,'MarkerEdgeColor','none','Color',gradientMap), %legend(string(rad2deg(TestedGrat_2)))
    plot(TestedDiff, abs(true-estim)); legend(string(norm_param_sigma))
    colororder(gradientMap)
    title('velocity error'),
    grid on,ylim([0,1.7])
end
% %% Debug function
% % How to Choose Parameters Systematically:
% % Step 1: Start with logistic_centre
% % 
% % Too high (>0.5): Kills most activity → single peak dominance
% % Too low (<0.1): Preserves too much noise
% % Sweet spot: 0.15 - 0.35 depending on your noise level
% % 
% % Step 2: Adjust K for competition strength
% % 
% % K = 1.5: Winner-take-all (only one peak survives)
% % K = 2-3: Moderate competition (2-3 peaks can coexist)
% % K > 3: Gentle competition (multiple peaks preserved)
% % 
% % Step 3: Tune spatial scales
% % 
% % Start with sigma_r and sigma_t roughly equal to your feature spacing
% % If your velocity grid spacing is 2 units, try sigma_r = 0.5 to sigma_r = 2
% % If orientation spacing is π/8, try sigma_t = π/16 to sigma_t = π/4
% % 
% % Step 4: Fine-tune logistic_slope
% % 
% % Low (1-3): Gentle contrast enhancement
% % Medium (4-7): Moderate contrast
% % High (8+): Aggressive thresholding
% 
% pop_resp = pop_resp_even;
% % If your data grid is, say, 8x11:
% sigma_r = (param(1).pref_vel(end) - param(1).pref_vel(end-1));
% sigma_t = pi/8;
% K = 2;
% logistic_slope = 5;
% logistic_centre = 0.4;
% max_iteration = 8;
% 
% debugDecodeMxHat(pop_resp(:,:,1,1,1),param(1),sigma_r,sigma_t,K,max_iteration,logistic_slope,logistic_centre);