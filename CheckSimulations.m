clear variables
close all

addpath FUNCTIONS\

folder2check = "NoNoise";

folder = dir(strcat("SIMULATIONS\PlaidAnalysis\New\",folder2check));

totFile = numel(folder) - 2;
weight_mode = 1;

for numFile = 1:totFile
    load(strcat("SIMULATIONS\PlaidAnalysis\New\",folder2check,"\",folder(numFile+2).name))
    for seeds = 1:25
    % M = 3e5; %compute from population response to plaid
    % M = max(data,[],"all");
    data = cat(5,totsim{1,:});
    stim = squeeze(cat(4,totsim{2,:}));
    param = squeeze(cat(4,totsim{3,:}));
    
    % num_or_ch_pooledTot = unique(cat(1,param.num_or_ch_pooled));
    num_diff_c = unique(cat(1,param.alpha));
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
    TotTestedDiff = (cat(1,stim.alpha)); 
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

    % 
    % plotPlaidVectors(repmat(testedTrueThetaPlaid,numel(TestedGrat_2),1),...
    %             squeeze(TotTestedGrat(1,:,1,1,:)),...
    %             repmat(testedVplaid,numel(TestedGrat_2),1), ...
    %             squeeze(TotTestedVGrat(1,:,1,1,:))) 
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
        weightingfunctions_general(xx,tt), ...
        1,1,1,prod(sze_reshaped(4:5)),sze_reshaped(6),sze_reshaped(7)), ...
        [1,2,4,5,6,3]);  
    
    %for each orientation combination make a figure
    %I plot what happen if I have the grating contrast level at different ratio
    %(columns) with different normalisation levels (rows)
    %EVEN POPULATION
    rng(seeds)
    for ii = 1:prod(sze_reshaped(4:5))
        M = max(squeeze(squeeze(data_reshaped(1,:,:,ii,:,:))),[],"all");
        n = (rand(size(squeeze(data_reshaped(1,:,:,ii,:,:)))))*M*.025;
        e_1(:,:,ii,:,:) = squeeze(data_reshaped(1,:,:,ii,:,:)) + n;
        M = max(squeeze(squeeze(data_reshaped(1,:,:,ii,:,:))),[],"all");
        n = (rand(size(squeeze(data_reshaped(2,:,:,ii,:,:)))))*M*.025;
        e_2(:,:,ii,:,:) = squeeze(data_reshaped(2,:,:,ii,:,:)) + n;
        
    
        w = squeeze(W(:,:,ii,:,:,:));
        pop_resp_odd(:,:,ii,:,:,:) = getpopresponse(w,e_1(:,:,ii,:,:));  % Apply weighting to the response data 
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
    sigma_r = (param(1).pref_vel(end) - param(1).pref_vel(end-1))*4/3;
    % sigma_r = 0.5;
    sigma_t = (pi/8)*0.8;
    K = 4;
    logistic_slope = 2;
    logistic_centre = 0.2;
    max_iteration = 2;
    
    pop_resp_even = squeeze(pop_resp_even(:,:,:,:,:,weight_mode));
    sze_pop = size(pop_resp_even);
    num_cond = prod(sze_pop(3:end));
    %Activity Decoding
    for ii = 1:num_cond
        [PR_decoded(:,:,:,ii),v_CM(:,ii),v_max(:,ii)] = DecodeMxHat( ...
                                    (squeeze(pop_resp_even(:,:,ii))) , ...  %Activity
                                    param(1), ...                                                              %Population Parameters
                                    sigma_r, ...                                                               %see code
                                    sigma_t, ...
                                    K, ...
                                    max_iteration, ...
                                    logistic_slope, ...
                                    logistic_centre);
    end
    % [M_decoded, ]  = max(reshape(pop_resp_even,[],num_cond),1);
   
    PR_decoded = reshape(PR_decoded,[sze_pop(1),sze_pop(2),max_iteration,sze_pop(3:end)]);
    
    %PLOT ACTIVITY AND DECODING
    if numFile == 1 && seeds == 1
        for numstim = 1:size(e_2,3)
            for normparam = 1:size(e_2,5)
                figure(50 + numstim + (normparam - 1)*size(e_2,3)), t = popresponse_tiled(cat(3,...
                                    reshape((squeeze(e_2(:,:,numstim,:,normparam))),[sze_pop(1:2),1,sze_pop(4)]),...
                                    reshape((squeeze(pop_resp_even(:,:,numstim,:,normparam))),[sze_pop(1:2),1,sze_pop(4)]),...
                                    reshape(squeeze(PR_decoded(:,:,:,numstim,:,normparam)),[sze_pop(1:2),max_iteration,sze_pop(4)])), ...
                                    param(1).pref_vel); 
                title(t,['PopActivity to stim ', num2str(numstim),' normalisation sigma ',num2str(norm_param_sigma(normparam))])
            end
        end
    end

    %Decode Activity with weighted activity
    for ii = 1:num_cond
        [vx(ii),vy(ii)] = decodingCM(squeeze(pop_resp_even(:,:,ii)),param(1));
    end
    % title(t,'Population Activity Decoded (rows are different Iteration, col are contrast)')
    vx = reshape(v_max(1,:),sze_pop(3:end));
    vy = reshape(v_max(2,:),sze_pop(3:end));
    gradientMap = lines(numel(norm_param_sigma));
    
    n_orient = 8;    
    x = linspace(1,n_orient,100);
    x_8 = round(linspace(1,100,n_orient));
    
    % figure(70),
    % w_o = exp(-(x(x_8(3))-(1:n_orient)).^2./(2.*norm_param_sigma.^2));
    % plot(w_o'),colororder(gradientMap), legend(string(norm_param_sigma)),
    % title('Tested Orientation Pooling (sample for one neuron)')
    % 
    v = sqrt(vx.^2 + vy.^2);
    ori = atan2(vy,vx);
    
    % FIX: Normalize angles to [0, 2π] for consistent comparison
    ori(ori < 0) = ori(ori < 0) + 2*pi;
    
    % Also ensure truetheta is in [0, 2π]
    true_theta_normalized = stim(1,1).truetheta;
    if true_theta_normalized < 0
        true_theta_normalized = true_theta_normalized + 2*pi;
    end
    
    for jj = 1:size(ori,1)
        figure(60),
        nexttile(jj),
        % subplot(2,1,1),
        hold on,
        estim = squeeze(ori(jj,:,:));
        true = repmat(true_theta_normalized, size(estim,1), size(estim,2));
        
        % Use angdiff for circular angular differences
        ang_error = abs(rad2deg(angdiff(estim, true)));
        
        scatter(TestedDiff+(rand(length(TestedDiff),1)-.5)*2*.007, ang_error, 15, 'filled', ...
            'Color', gradientMap,'LineWidth',0.3,'MarkerFaceAlpha',0.15,'MarkerEdgeColor','none');
        % p = plot(TestedDiff, ang_error, 'LineStyle','-.');
        % for pp = 1:numel(p)
        %     p(pp).Color(4) = 0.3;
        % end

        legend(string(norm_param_sigma))
        % colororder(gradientMap)
        ylim([0,100]),xlim([0,max(TestedDiff)]) 
        grid on
        title('angle error')
        ylabel("[deg]"), xlabel("[delta-contrast]")
        tot_ang_error(:,:,jj,seeds) = ang_error;
        % % velocity
        % estim = squeeze(v(jj,:,:,weight_mode));
        % true = repmat(stim(1,1).vpld, size(estim,1), size(estim,2));
        % subplot(2,1,2),
        % hold on,
        % scatter(TestedDiff, abs(true-estim), 'filled', ...
        %     'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'none', 'Color', gradientMap);
        % plot(TestedDiff, abs(true-estim)); 
        % legend(string(norm_param_sigma))
        % colororder(gradientMap)
        % title('velocity error'),
        % grid on, 
        % ylim([0,5])
        % ylabel("[pix/frames]"), xlabel("[delta-contrast]")
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
    clear pop_resp_even pop_resp_odd pop_resp
    end
    figure(60)
     for jj = 1:size(ori,1)
        figure(60),
        nexttile(jj),
        % subplot(2,1,1),
        hold on,
        errorbar(TestedDiff, squeeze(mean(tot_ang_error(:,:,jj,:),4)),squeeze(std(tot_ang_error(:,:,jj,:),[],4)))
        legend(string(norm_param_sigma))
        colororder(gradientMap)
        % colororder(gradientMap)
     end
end

% v_est = reshape(v_max,[2,sze_pop(3:end)]);
% % close all
% 
% %visualize estimates
% for weight_mode = 2
% for norm = 1
% figure, 
% for numstim = 1:3
% nexttile(numstim),
% scatter(squeeze(v_est(1,numstim,:,norm,weight_mode)),squeeze(v_est(2,numstim,:,norm,weight_mode)),'filled')
% xlim([-4,4]), ylim([-4,4]), 
% grid on
% end
% end
% end