clear variables
load("SIMULATIONS\PlaidAnalysis\SimulationTot_NormEffect20250616_180425.mat")

data = cat(5,totsim{1,:});
sze = size(data);
numComplexCellsPopulations = sze(1);
numOrient = sze(2);
numVel = sze(3);
numTestedStim = sze(4);       %stim parameters tested
numTestedParameters = sze(5); %population parameters
stim = squeeze(cat(4,totsim{2,:}));
param = squeeze(cat(4,totsim{3,:}));

%stim tested parameters
testedTheta_grat = rad2deg(permute(cat(3,stim.theta_g),[1,3,2]));

testedContrasts = permute(cat(3,stim.contrast),[1,3,2]);
testedTruethetaPlaid = rad2deg(cat(2,stim.truetheta));
%population tested parameters
norm_param = cat(1,param.norm_param);

%Look at data
figure,stim(1).disp = 1;
pl = plaid(stim(1));
generate_plaid(pl);
figure, plotPopResponse(squeeze(data(1,:,:,1,1)),param(1).pref_vel)
title(['Plaid at \theta = ',num2str(testedTruethetaPlaid(1,1)), ...
    ' grat orient = ','[',num2str(testedTheta_grat(1,1,1)),',',num2str(testedTheta_grat(1,1,2)),']', ...
    sprintf('\nContrast diff = '), ...
    '[',num2str(testedContrasts(1,1,1)),',',num2str(testedContrasts(1,1,2)),']'])

