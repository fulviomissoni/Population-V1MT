function C1 = popFlowV1MT(II,param)
% This function computes a distributed analysis of a stack of visual images 
% by means of a population of motion energy detectors. 
%
% -) [rows,cols,n_frames] = size(II); represents the stereo-input
% -) 'thr' contains the threshold values below which activity
%       at rest is assumed -> in particular are related to the two channels of
%       representation: "energy" and "orientation"
% -) parameters is a 7x1 struct and contains all the specifics of the filters
%     param.spat_freq: spatial frequency of spatial Gabor filter
%     param.n_orient:  number of orientation channels
%     param.pref_vel:  vector of preferred velocity in [pix/frame]
%     param.temp_filt: model of temporal component of the spatio-temporal RFs ('gabor','exp_decay','adelson_bergen')
%     param.spatial_filt: name (in char) of the file that contains N_o 1D gabor filters that are necessary to construct the 2D oriented Gabor filters
%     param.samples: number of samples of the spatial filter
%     param.norm_param: normalization factors used in the two stages (layer 2 and layer 3)
%     param.sigma_pool: sigma values of gaussian for spatial pooling;
%     param.num_or_ch_pooled: number of orientation channels pooled; 

%%%%%%%%%%%%%%

k0 = param.spat_freq;
v = param.pref_vel;
Ft_choice = param.temp_filt;
filter_file = param.spatial_filt;

% COMPUTE GABOR FILTERING (SPATIAL AND TEMPORAL PART)

if param.n_orient == 8
    [F] = filtGaborSpace2D(II,filter_file);
end
if param.n_orient == 16
    [F] = filtSepGabor(II,filter_file);
end

F = filtTime(F,'valid',Ft_choice,v,k0);

for i=1:4
    %reverse position of third and fourth dimension
    F{i} = permute(F{i},[1 2 4 3 5]);
end

%[nr, nc, n_frame, or, dumb] = size(F{1});

% POPULATION ENCODING
% This function combines oppurtunely filtered responses to obtain
% population activity of V1-like cells and motion energy detectors (C1)
C1 = populationV1MT(F,param);

clear F


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C1] = populationV1MT(G,param)

% This function combines oppurtunely filtered responses to obtain
% population activity of V1-like cells and motion energy detectors (C1)
alpha = param.norm_param;
sigma_pool = param.sigma_pool;
num_or_ch_pooled = param.num_or_ch_pooled;
w_o = param.orient_weighting;

[sy, sx, n_frames, n_orient, v] = size(G{1});
sze = size(G{1});

C = G{1}(:);    % REAL LEFT
S = G{2}(:);    % IMAG LEFT    
Ct = G{3}(:);   % REAL LEFT - Temporal Derivative
St = G{4}(:);   % IMAG LEFT - Temporal Derivative

clear G;
w_orient = permute(repmat(w_o,sy*sx*n_frames,1,v),[2,3,1]);
%% ENERGY MODEL

% ALLOCATE MEMORY
S0 = cell(8,1);
for i=1:4
    S0{i} = zeros(sy*sx*n_frames*n_orient*v,1);
end

%SPATIO TEMPORAL ORIENTED FILTERS
S0{1} = S+Ct;
S0{2} = -St+C;
S0{3} = -S+Ct;
S0{4} = St+C;
clear C S Ct St

% COMPLEX CELLS - FIRST LAYER
C1{1} = S0{1}.^2 + S0{2}.^2;
C1{2} = S0{3}.^2 + S0{4}.^2;

clear S0
% THRESHOLDING
th_C(1) = median(C1{1}(:));th_C(2) = median(C1{2}(:));

C1{1}(C1{1}(:)<0.01*th_C(1)) = 0;
C1{2}(C1{2}(:)<0.01*th_C(2)) = 0;

% NORMALIZATION STAGE OF COMPLEX-CELLS
a1 = alpha(1);
a2 = alpha(2);

%Memory allocation:

for i = 1:2
    C1_pooled{i} = reshape(C1{i},sy,sx,n_frames*n_orient*v);
    S = zeros(sy,sx,n_frames*n_orient*v);
    %spatial pooling of normalization pool
    for p = 1:n_frames*n_orient*v
        tmp = C1_pooled{i}(:,:,p);
        % tmp2 = imgaussfilt(tmp,sigma_pool);
        tmp2 = tmp;
        S(:,:,p) = tmp2;
        C1_pooled{i}(:,:,p) = tmp2;
    end
    % C1_pooled{i} = reshape(C1_pooled{i},sy*sx*n_frames,n_orient,v);
    C1_pooled{i} = reshape(C1_pooled{i},1,[]);

    S = reshape(S,sy*sx*n_frames,n_orient,v);
    S = permute(S,[2 3 1]);
    % m = 1/mean(S,'all');
    m = 1;
    % eps = 0.025*mean(S,'all');
    %orientation pooling
    S = repmat(sum(pagemtimes(w_o, S),2),[1 v 1]); %apply orientation weighting
    
%     a2 = 1;
%     a2 = 1/max(max(S));

    %orientation pooling

%     % if num_or_ch_pooled == 8
%     %     S = sum(S.*w_orient,[1,2]);
%     %     S = repmat(S,[n_orient v 1]);
%     % end
%     % if num_or_ch_pooled == 1
%     %     S = sum(S,2);
%     %     S = repmat(S,[1 v 1]);
%     % end
% %     if num_or_ch_pooled~=8
% %         index_o = circshift(index_o,floor(num_or_ch_pooled/2));
% %         for o = 1:n_orient  
% %             S(o,:) = sum(tmp(index_o(1:1+num_or_ch_pooled-1),:));
% %             index_o = circshift(index_o,-floor(num_or_ch_pooled/2));
% %         end
% %     end
% %     if num_or_ch_pooled==8
% %         S = repmat(sum(tmp),[n_orient,1]);
% %     end
% %     if num_or_ch_pooled==1
% %         S = tmp;
% %     end
    % S = reshape(S,n_orient,v,sy*sx*n_frames);
    S = permute(S,[3 1 2]);
    S = reshape(S,1,[]);
    %THRESHOLDING S
    S(S<0.3*median(S(:))) = 1e-19;
    C1_pooled{i}(S<0.3*median(S(:))) = 1e-19;
    % a1 = 0.001*mean(C1_pooled{i},'all');
    % a2 = a1 / mean(S(:));
    a1 = 1; a2 = 1;
    C1{i} = C1_pooled{i}./(a1 + a2*S);
    % orientation_activity = repmat(squeeze( ...
    %     prctile( ...
    %     reshape(C1_pooled{i},sze), ...
    %     90,[1,2,5]))',[sze(1),1,sze(2),sze(3),sze(5)]);
    % orientation_activity  = permute(orientation_activity,[1 , 3, 4, 2, 5]);
    % global_activity = prctile(C1_pooled{i}(:),95);
    % C1{1} = C1_pooled{i}./(a1+ + 0.005*global_activity + S);
    % C1{i} = C1_pooled{i}./(orientation_activity(:)' + S);
    C1{i}(isnan(C1{i})) = 0;
    % C1{i}(C1)
end

for i=1:2
    C1{i} = squeeze(reshape(C1{i},sze));
end