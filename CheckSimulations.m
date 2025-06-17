clear variables
close all
load("SIMULATIONS\PlaidAnalysis\NoNoise\SimulationTot_NormEffect20250617_160441.mat")
% load("SIMULATIONS\PlaidAnalysis\Noise\SimulationTot_Noise_NormEffect20250617_160441.mat")

% M = 3e5; %compute from population response to plaid
% M = max(data,[],"all");
data = cat(5,totsim{1,:});
stim = squeeze(cat(4,totsim{2,:}));
param = squeeze(cat(4,totsim{3,:}));

num_or_ch_pooledTot = unique(cat(1,param.num_or_ch_pooled));
num_diff_c = unique(cat(1,param.diff_c));

sze = size(data);
numComplexCellsPopulations = sze(1);
numOrient = sze(2);
numVel = sze(3);
numTestedStim = sze(4);       %stim parameters tested
numTestedParameters = sze(5); %population parameters



data_reshaped = reshape(data,[sze(1:end-1),numel(num_or_ch_pooledTot),numel(num_diff_c)]);
param_reshaped = reshape(param,[numel(num_or_ch_pooledTot),numel(num_diff_c)]);

% data = logistic(cat(5,totsim{1,:}),2.5e2,.003);


%stim tested parameters
testedTheta_grat = rad2deg( ...
    reshape( ...
    cat(1,stim.theta_g), ...
    [size(stim,1),numel(num_or_ch_pooledTot),numel(num_diff_c),2]));
testedV_grat = abs( ...
    reshape( ...
    cat(1,stim.vgrat), ...
    [size(stim,1),numel(num_or_ch_pooledTot),numel(num_diff_c),2]));
testedContrasts = reshape( ...
    cat(1,stim.contrast), ...
    [size(stim,1),numel(num_or_ch_pooledTot),numel(num_diff_c),2]);
testedTrueThetaPlaid = rad2deg( ...
    reshape( ...
    cat(1,stim.truetheta), ...
    [size(stim,1),numel(num_or_ch_pooledTot),numel(num_diff_c)]));
%population tested parameters
norm_param = cat(1,param.norm_param);
plotTestedVelocitiesDirection

%% Look at the data

% figure,
% stim(1).disp = 1;
% pl = plaid(stim(1));
% generate_plaid(pl);
for ii = 1:numel(num_or_ch_pooledTot)
    figure,
    for jj = 1:numel(num_diff_c)
        % j = mod(ii-1, numel(param.pref_vel)) + 1;
        subplot(1,numel(num_diff_c),jj),

        dataplot = squeeze(data_reshaped(2,:,:,1,ii,jj));
        % Med = median(dataplot,"all");
        % dataplot = logistic(dataplot,mean(dataplot),.002);
        plotPopResponse(dataplot,param_reshaped(ii,jj).pref_vel)
        title(['Plaid at \theta = ',num2str(testedTrueThetaPlaid(1,ii,jj)), ...
            ' grat orient = ','[',num2str(testedTheta_grat(1,ii,jj,1)),',',num2str(testedTheta_grat(1,ii,jj,2)),']', ...
            sprintf('\nContrast diff = '), ...
            '[',num2str(testedContrasts(1,ii,jj,1)),',',num2str(testedContrasts(1,ii,jj,2)),']'])
        % set(gca,'CLim',[0,1.1])
        
        %% Test Weighting 
        
        theta_cell_OUT = 0:pi/param_reshaped(ii,jj).n_orient:pi-pi/param_reshaped(ii,jj).n_orient;    
        [xx,tt] = meshgrid(param_reshaped(ii,jj).pref_vel,theta_cell_OUT);
        sigma = 0.25;
        W2 = exp( ...
            -(xx(:).*cos(tt(:)-tt(:)') - xx(:)').^2./ ...
            (2*sigma.^2)); %Andre pesi originali (inviluppo exp(.))
        % W2 = *cos(tt(:)-tt(:)'); %Chessa-Solari model
        % W2 = cos(xx(:).*cos(tt(:)'-tt(:)) - xx(:)'); %primo modello per introdurre l'opponenza (inviluppo coseno)
        % W2 = exp(-(xx(:).*cos(tt(:)'-tt(:)) - xx(:)').^2./(2*sigmaGabor.^2)).* ...
        % cos(2*pi*1./(4*xx(:)').*(xx(:).*cos(tt(:)'-tt(:)) - xx(:)')); %inviluppo Gabor per gestire la quantità e la posizione dei pesi negativi
        % W2 = reshape(reshape(W2,8,11,8,11)./max(reshape(W2,8,11,8,11),[],4),88,88);
        W2(isnan(W2)) = 0;
        % pop_resp = reshape(pop_resp,param.n_orient,numel(param.pref_vel),numel(param.num_or_ch_pooled),numel(param.diff_c));
        
        dataplot(dataplot<0) = 0;
        PR_norm = dataplot./mean(dataplot,"all");
        
        MT_norm = reshape(W2*reshape(PR_norm,sze(2)*sze(3),[]),sze(2:3));
        % nexttile(jj+numel(num_diff_c))
        % plotPopResponse(MT_norm,param_reshaped(1).pref_vel)
    end
end