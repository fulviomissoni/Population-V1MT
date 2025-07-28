%% %METODO FLIP
close all

% full orientation pool
%P is activity response with spatial pooling
%consider 16 orientation channels
PFull = [P(:,6:end,:);fliplr(P(:,1:6,:))];
%pool all channels
PFull = repmat(sum(PFull),[16,1,1]);
%re-organized matrices
% Pool = [PFull(1:8,2:end,:),fliplr(PFull(9:end,1:end,:))];
Pool = [fliplr(PFull(9:16,2:end)),PFull(1:8,:)];
%normalize
CTNorm =CT./(1+1/(mean(mean(Pool))).*Pool);
%plot data
figure, plotPopResponse(squeeze(CT(:,:,1)),0,0,param.pref_vel)
figure, plotPopResponse(squeeze(CTNorm(:,:,1)),0,0,param.pref_vel)
figure, plotPopResponse(squeeze(Pool(:,:,1)),0,0,param.pref_vel)
%%
% single orientation channels
%consider 16 orientation channels
Pool = CT;
%normalize
CTNorm = CT./(1+1/(mean(mean(Pool)))*Pool);
%plot data
figure, plotPopResponse(CT,0,0,param.pref_vel)
figure, plotPopResponse(CTNorm,0,0,param.pref_vel)
figure, plotPopResponse(Pool,0,0,param.pref_vel)
%% SUM along velocities
% all orientations 
P = CT;
PFull = sum(sum(P));
Pool = PFull;
%normalize
CTNorm =CT./(1+1/(mean(mean(Pool)))*Pool);
%plot data
figure, plotPopResponse(CT,0,0,param.pref_vel)
figure, plotPopResponse(CTNorm,0,0,param.pref_vel)
figure, plotPopResponse(repmat(Pool,[8,11]),0,0,param.pref_vel)
% single orientation
P = CT;
PFull = sum(P,2);
Pool = repmat(PFull,[1,11]);
%normalize
CTNorm =CT./(1+1/(mean(mean(Pool)))*Pool);
%plot data
figure, plotPopResponse(CT,0,0,param.pref_vel)
figure, plotPopResponse(CTNorm,0,0,param.pref_vel)
figure, plotPopResponse(Pool,0,0,param.pref_vel)