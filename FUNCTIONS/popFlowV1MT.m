function C1 = popFlowV1MT(II,param)
% This function computes a distributed analysis of a stack of visual images 
% by means of a population of motion energy detectors. Th
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

[sy, sx, n_frames, n_orient, v] = size(G{1});
sze = size(G{1});

C = G{1}(:);    % REAL LEFT
S = G{2}(:);    % IMAG LEFT    
Ct = G{3}(:);   % REAL LEFT - Temporal Derivative
St = G{4}(:);   % IMAG LEFT - Temporal Derivative

clear G;

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

% NORMALIZATION STAGE OF COMPLEX-CELLS
a1 = alpha(1);
a2 = alpha(2);

for i = 1:2
    C1{i} = reshape(C1{i},sy,sx,n_frames*n_orient*v);
    S = zeros(sy,sx,n_frames*n_orient*v);
    %spatial pooling of normalization pool
    for p = 1:n_frames*n_orient*v
        tmp = C1{i}(:,:,p);
        tmp2 = imgaussfilt(tmp,sigma_pool);
        S(:,:,p) = tmp2;
    end
    C1{i} = reshape(C1{i},1,[]);
    S = reshape(S,sy*sx*n_frames,n_orient,v);
    S = permute(S,[2 3 1]);
%     S = reshape(S,n_orient,[]);
    index_o = 1:8;
    tmp = S;
    m = 1/mean(mean(mean((S))));
%     a2 = 1;
%     a2 = 1/max(max(S));

    %orientation pooling
    if num_or_ch_pooled == 8
        S = sum(sum(S,2),1);
        S = repmat(S,[n_orient v 1]);
    end
    if num_or_ch_pooled == 1
        S = sum(S,2);
        S = repmat(S,[1 v 1]);
    end
%     if num_or_ch_pooled~=8
%         index_o = circshift(index_o,floor(num_or_ch_pooled/2));
%         for o = 1:n_orient  
%             S(o,:) = sum(tmp(index_o(1:1+num_or_ch_pooled-1),:));
%             index_o = circshift(index_o,-floor(num_or_ch_pooled/2));
%         end
%     end
%     if num_or_ch_pooled==8
%         S = repmat(sum(tmp),[n_orient,1]);
%     end
%     if num_or_ch_pooled==1
%         S = tmp;
%     end
    S = reshape(S,n_orient,v,sy*sx*n_frames);
    S = permute(S,[3 1 2]);
    S = reshape(S,1,[]);
    
    C1{i} = C1{i}./(a1 + a2*m*S);
    C1{i}(isnan(C1{i})) = 0;
end

for i=1:2
    C1{i} = squeeze(reshape(C1{i},sze));
end