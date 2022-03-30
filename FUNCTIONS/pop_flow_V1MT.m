function [MT, EC21, EC22, varargout] = pop_flow_V1MT(II,param,tun_flag)
% This function computes a distributed analysis of a stereo video input with two population of
% descriptors organized hierarchically. These detectors are constructed as
% described in my thesis. The first population consists in a set of
% disparity detectors that provides distributed representation of the
% static disparity map. This representation is coded by the population
% activity. The second population, instead, consists in a set of more
% complex detectors that are selective for the temporal variation of the
% disparity.
%
% -) [rows,cols,n_frames,2] = size(II); represents the stereo-input
% -) 'thr' contains the threshold values below which activity
%       at rest is assumed -> in particular are related to the two channels of
%       representation: "energy" and "orientation"
% -) parameters is a 7x1 struct and contains all the specifics of the filters
%     param.spatFreq: spatial frequency of spatial Gabor filter
%     param.ocDom:    ocular-dominance
%     param.phShift:  phase_shifts for each orientation channel
%     param.nOrient:  number of orientation channels
%     param.prefVel:  vector of preferred velocity in [pix/frame]
%     param.tempFilt: model of temporal component of the spatio-temporal RFs ('gabor','exp_decay','adelson_bergen')
%     param.spatialFilt: name (in char) of the file that contains N_o 1D gabor filters that are necessary to construct the 2D oriented Gabor filters
%     param.sample: number of samples of the spatial filter
%     param.normParam: normalization factors used in the two stages (layer 2 and layer 3)
%     

%%%%%%%%%%%%%%

k0 = param.spatFreq;
a = param.ocDom;
ph = param.phShift;
v = param.prefVel;
Ft_choice = param.tempFilt;
filter_file = param.spatialFilt;
alpha = param.normParam;


% COMPUTE GABOR FILTERING (SPATIAL AND TEMPORAL PART)
% tic
if param.nOrient == 8
    [F] = filt_gabor_space2D(II,filter_file,tun_flag);
end
if param.nOrient == 16
    [F] = filt_sep_gabor(II,filter_file);
end
% toc
%if tun == 1 then consider only the cell in the center of the frame
if tun_flag == 1 
    [sy,sx,dumb] = size(F{1,1});
    F{1} = F{1}(ceil(sy/2)-1,ceil(sx/2)-1,:,:,:,:);
    F{2} = F{2}(ceil(sy/2)-1,ceil(sx/2)-1,:,:,:,:);
end
% tic
F = filt_time(F,'valid',Ft_choice,v,k0);
% toc
for i=1:4
    %move side channel on to the sixth dimension
    F{i} = permute(F{i},[1 2 4 3 6 5]);
end
[nr, nc, n_frame, or, dumb] = size(F{1});
% ALLOCATE MEMORY
G = cell(4,2);
for i=1:4
    G{i,1} = F{i}(:,:,:,:,:,1);
    G{i,2} = zeros(nr,nc,n_frame,or,length(v),size(ph,2));
end
% COMPUTE PHASE SHIFT 
[G{1,2},G{2,2}] = shift_in_phase(F{1}(:,:,:,:,:,2),F{2}(:,:,:,:,:,2),ph);  
[G{3,2},G{4,2}] = shift_in_phase(F{3}(:,:,:,:,:,2),F{4}(:,:,:,:,:,2),ph); 
clear  F
%G contains the distributed response of a population of monocular separable 
%recptive fields 

% POPULATION ENCODING
% This function combine opportunely these responses to obtain the disparity
% detectors response and the MT-like population response
[MT, EC21, EC22, C1,S] = populationV1MT(G,a,alpha);
MT(isnan(MT)) = 0;
EC21(isnan(EC21)) = 0;
EC22(isnan(EC22)) = 0;
varargout{1} = C1;
varargout{2} = S;

clear G


function [MT, EC21, EC22, C1,varargout] = populationV1MT(G,a,varargin)


[sy, sx, n_frames, n_orient, v, phase_num] = size(G{1,2});

CL = G{1,1}(:);    % REAL LEFT
SL = G{2,1}(:);    % IMAG LEFT    
CtL = G{3,1}(:);   % REAL LEFT - Temporal Derivative
StL = G{4,1}(:);   % IMAG LEFT - Temporal Derivative

CR = G{1,2};    % REAL RIGHT
CR = reshape(CR,[],phase_num);
SR = G{2,2};    % IMAG RIGHT
SR = reshape(SR,[],phase_num);
CtR = G{3,2};   % REAL RIGHT - Temporal Derivative
CtR = reshape(CtR,[],phase_num);
StR = G{4,2};   % IMAG RIGHT - Temporal Derivative
StR = reshape(StR,[],phase_num);

clear G;
%% ENERGY MODEL
% ALLOCATE MEMORY FOR ENERGY
S0 = cell(8,1);
for i=1:4
    S0{i} = zeros(sy*sx*n_frames*n_orient*v,1);
    S0{i+4} = zeros(sy*sx*n_frames*n_orient*v,phase_num);
end
% LEFT CHANNEL
S0{1} = SL+CtL;
S0{2} = -StL+CL;
S0{3} = -SL+CtL;
S0{4} = StL+CL;
clear CL SL CtL StL
% RIGHT CHANNEL
S0{5} = -StR+CR;
S0{6} = SR+CtR;
S0{7} = StR+CR;
S0{8} = -SR+CtR;
clear CR SR CtR StR

% COMPLEX CELLS - FIRST LAYER
for ph_n=1:phase_num
    S1{1}(:,ph_n) = ((1-a).*S0{1}-a.*squeeze(S0{5}(:,1)));
    S1{2}(:,ph_n) = ((1-a).*S0{2}+a.*squeeze(S0{6}(:,1)));
    S1{3}(:,ph_n) = ((1-a).*S0{3}-a.*squeeze(S0{7}(:,1)));
    S1{4}(:,ph_n) = ((1-a).*S0{4}+a.*squeeze(S0{8}(:,1)));
    C1{1}(:,ph_n) = S1{1}.^2+S1{2}.^2;
    C1{2}(:,ph_n) = S1{3}.^2+S1{4}.^2;
    C1{3}(:,ph_n) = (a.*S0{1}-(1-a).*squeeze(S0{5}(:,1))).^2+(a.*S0{2}+(1-a).*squeeze(S0{6}(:,1))).^2;
    S0{5}(:,1) = [];
    S0{6}(:,1) = [];
    C1{4}(:,ph_n) = (a.*S0{3}-(1-a).*squeeze(S0{7}(:,1))).^2+(a.*S0{4}+(1-a).*squeeze(S0{8}(:,1))).^2;
    S0{7}(:,1) = [];
    S0{8}(:,1) = [];
end

clear S0
% NORMALIZATION STAGE OF COMPLEX-CELLS
if nargin>2
    a1 = varargin{1}(1,1);
    a2 = varargin{1}(2,1);
else
    a1 = 1; a2 = 0;
end

sze = size(C1{1});
for i = 1:4
    C1{i} = reshape(C1{i},sy,sx,n_frames*n_orient*v*phase_num);
    S = zeros(sy,sx,n_frames*n_orient*v*phase_num);
%     S = squeeze(mean(mean(C1{i})));
    sigmaPool = 3;
    %spatial pooling of normalization pool
    for p = 1:n_frames*n_orient*v*phase_num
        tmp = C1{i}(:,:,p);
        tmp2 = imgaussfilt(tmp,sigmaPool);
        S(:,:,p) = tmp2;
    end
    C1{i} = reshape(C1{i},sy*sx*n_frames*n_orient,v*phase_num);
    S = reshape(S,sy*sx*n_frames,n_orient,v*phase_num);
    S = permute(S,[2 1 3]);
    S = reshape(S,n_orient,[]);
    index_o = circshift(1:n_orient,3);
    %orientation pooling
    for o = 1:n_orient
        tmp = S(index_o(1:5),:);
        S(o,:) = sum(tmp);
        index_o = circshift(index_o,-1);
    end
    S = reshape(S,n_orient,sy*sx*n_frames,v*phase_num);
    S = permute(S,[2 1 3]);
    S = reshape(C1{i},sy*sx*n_frames*n_orient,v*phase_num);
    C1{i} = C1{i}./(a1 + a2*S);
    C1{i}(isnan(C1{i})) = 0;
%     C1{i} = reshape(C1{i},sze);
end

%connection weights
% copp = 0.6;
copp = 1;
% COMPLEX CELLS - SECOND LAYER
C21 = C1{2}-copp*C1{1}; C22=C1{3}-copp*C1{4};
% % Half-wave rectification of the response (is a Firing Rate model)
% index = C21<0;
% C21(index) = 0;
% index = C22<0;
% C22(index) = 0;
%reshape for the output
for i=1:4
    C1{i} = reshape(C1{i},sy,sx,n_orient,n_frames,v,phase_num);
end
%Disparity Detectors response OUTPUT
EC21 = reshape(C21,sy,sx,n_orient,n_frames,v,phase_num);
EC22 = reshape(C22,sy,sx,n_orient,n_frames,v,phase_num);
EC21(isnan(EC21)) = 0;
EC22(isnan(EC22)) = 0;

% COMPLEX CELLS - LAYER THREE
C3 = C21 + C22;

% % NORMALIZATION of C3-cells
% if nargin>2
%     a1 = varargin{1}(1,2);
%     a2 = varargin{1}(2,2);
% else
%     a1 = 0; a2 = 1;
% end

% %permute to reduce the numbers of nested for and, thus, the computational time
% C3 = reshape(C3,sy*sx,n_orient,n_frames*v,phase_num);
% C3 = permute(C3,[1 3 2 4]);
% C3 = reshape(C3,sy*sx*n_frames*v,n_orient*phase_num);
% S = zeros(sy*sx*n_frames*v,1);
% for maps=1:n_orient*phase_num
%     S = S + C3(:,maps);
% end
% S = S/(n_orient); 
% for maps=1:n_orient*phase_num
%     C3(:,maps) = C3(:,maps)./(a1*S + a2);
% end
% C3 = reshape(C3,sy,sx,n_frames,v,n_orient,phase_num);
% C3 = permute(C3,[1 2 5 3 4 6]);

%%MT LAYER

% SPATIAL POOLING
C3 = reshape(C3,sy,sx,n_frames*n_orient*v*phase_num);
%Gaussian filter
filt = fspecial('gaussian', [5 5], 1);
% filt=filt./sum(sum(fspecial('gaussian', [31 31], 1))); %it's already
% normalized
for maps=1:(n_frames*n_orient*phase_num)
    C3(:,:,maps) = conv2b(C3(:,:,maps),filt,3);
end
C3 = reshape(C3,sy,sx,n_frames,n_orient,[]);
C3 = permute(C3,[1 2 4 3 5]);
C3 = reshape(C3,sy,sx,n_orient,n_frames*v*phase_num);
clear filt

% ORIENTATION POOLING
csum = zeros(sy,sx);
MT = zeros(sy,sx,n_orient,n_frames*v*phase_num);
for maps=1:n_frames*v*phase_num
    for d=1:n_orient
        thetad = (d-1)*pi/8; %the direction channel are at same of the orientation channel
        for o=1:n_orient
            theta = (o-1)*pi/8;
            csum = csum + C3(:,:,o,maps)*cos(theta-thetad);  
        end
        MT(:,:,d,maps) = csum;
        csum = zeros(sy,sx);
    end
end

MT = reshape(MT,sy,sx,n_orient,n_frames,v,phase_num);
clear C3
% %remove element<0 
% index = MT<0;
% MT(index) = 0;
MT(isnan(MT)) = 0;
varargout{1} = reshape(cat(1,S1{1},S1{2},S1{3},S1{4}),4,sy,sx,n_orient,n_frames,v,phase_num);