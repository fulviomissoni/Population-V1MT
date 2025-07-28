clear 
close all
clc

load('SIMULATIONS/PlaidAnalysis/vel_tuning_All_PlaidII_lambda0.mat')
load('SIMULATIONS/RDSAnalysis/vel_tuning_polarRDS_dur72.mat')
load('SIMULATIONS/PlaidAnalysis/pop_resp.mat')
load('SIMULATIONS/MTWeights.mat')

% %select one cell and the corresponding combination of gratings from EC
% %(All PlaidII data)
% TC = squeeze(EC(3,:,:,1,1,:,:)); %TUNING CURVES OF ALL V1/MT COMPONENT CELLS
TC = squeeze(e(3,:,:,:,:));
% TC = permute(pop_resp,[2 3 1]);
% TC = reshape(TC,param.n_orient,numel(param.pref_vel),numel(param.num_or_ch_pooled),numel(param.diff_c));
sze = size(TC);
%controlliamo di aver scelto giusto (non ho un valore molto simile ma
%lasciamo la coppia (1,1)
figure, plotPopResponse(squeeze(TC(:,:,3,10)),param.pref_vel)
%create MT pattern cells tuning curves
theta_cell_OUT = 0:pi/param.n_orient:pi-pi/param.n_orient;    
[xx,tt] = meshgrid(param.pref_vel,theta_cell_OUT);
sigma = 2;
W2 = exp(-(xx(:).*cos(tt(:)'-tt(:)) - xx(:)').^2/(2*sigma^2)).*cos(xx(:).*cos(tt(:)'-tt(:)) - xx(:)');
W2 = W2 - eye(size(W2));
W = reshape(W,sze(1)*sze(2),sze(1)*sze(2));
MT_TC = W*reshape(TC,sze(1)*sze(2),[]);
% MOVSHON
% Devo costruirmi la likelihood della popolazione in risposta ad un plaid 
% di tipo II che mi interessa - la costruisco sommando le
% log delle curve di tuning pesate dall'attività di popolazione a questo
% plaid (ottengo la log likelihood in realtà)
for i=1:numel(param.num_or_ch_pooled)
    for j=1:numel(param.diff_c)
        %select one stimulus
        pop_resp_BioGautama = reshape(pop_resp_BioGautama,numel(param.num_or_ch_pooled),numel(param.diff_c),param.n_orient,numel(param.pref_vel));
        PRMT = squeeze(pop_resp_BioGautama(i,j,:,:));
        LL = PRMT(:)'*log(reshape(TC,sze(3)*sze(4),[])');
        LL(LL<0) = 0;
        figure, plotPopResponse(reshape(LL,param.n_orient,numel(param.pref_vel)),param.pref_vel)
        pause
    end
end