function [varargout] = motion_popV1MT(param,stim)
% % % % % % % DEBUG % % % % % % %
debug = false;

file_name='mono_motion_pop';
samples=param.samples;
k0=param.spatFreq;
dur = stim.dur;
if ~isfield(stim,'type')
    error('Define type of the stimulus!! Options are: grat, plaid, or RDS')
end
fieldNames{1} = 'dur'; 
fieldNames{2} = 'vgrat'; fieldNames{3} = 'theta_g';
% fieldNames{4} = 'truetheta'; fieldNames{5} = 'contrast_g';
% fieldNames{6} = 'mode'; fieldNames{7} = 'vel_stim';
argCheck(stim,fieldNames);
%% Stimulus definition 
stimuli = ["plaid","grat","RDS_tuning","RDS_moving","shift_grat"];
if ~sum(matches(stimuli,stim.type))
    error('Select from the possible stimuli:\n %s, %s, %s',stimuli{1},stimuli{2},stimuli{3})
end
    
switch stim.type
    case 'grat'
        II = sinGrating(samples,samples,dur,[0,0],stim.vgrat,k0,stim.theta_g);
    case 'plaid'
        if ~isfield(stim,'pl_type')
            pl_type = 1;
        else 
            pl_type = stim.pl_type;
        end
        
        for num_pld = 1:length(stim.truetheta)
            %plaid object
            arg.dur = dur;                              %aperture size      [pixs]
            arg.apert_rad = ceil(samples/2)+2;          %aperture size      [pixs]
            arg.truetheta = stim.truetheta(num_pld);    %true orientation   [rad]
            arg.vpld = stim.vel_stim(num_pld);          %velocity amplitude [pixs/frame]
            arg.k = [k0,k0];                            %spatial freq       [cycle/pix]
            arg.vgrat = stim.vgrat(num_pld,:);          %gratings vel       [pixs/frame]
            arg.theta_g = stim.theta_g(num_pld,:);      %gratings orient    [rad]
            arg.alpha = 0.5;                            %alpha channel for transparency
            arg.contrast = stim.contrast_g(num_pld,:);             %Contrast of two gratings
            arg.mode = stim.mode;                       %stimulus implementation algorithm
            arg.pl_type = pl_type;                      %plaid type
            %define plaid object
            II{num_pld} = plaid(arg);
            %generate plaid stimulus
        end
    case 'RDS_tuning'         
        for num_stim = 1:length(stim.truetheta)
                vx = stim.vgrat(num_stim,1);
                vy = stim.vgrat(num_stim,2);
                II{num_stim} = moving_RDS_MS(samples,samples,dur,4,vx, vy);
                II{num_stim} = II{num_stim}(60:end-60,60:end-60,:);
        end
        II = reshape(II,stim.stim_size);
    case 'RDS_moving'   
        II{1} = moving_RDS_MS(samples,samples,43,4,stim.vgrat(1),stim.vgrat(2));
        
        %make moving dots
    case 'shift_grat'
        II = sinGrating(samples,samples,dur,[0,0],stim.vgrat(1),k0,stim.theta_g);
        for i=1:floor(1/(2*k0))+1
           II(i+1) = sinGrating(samples,samples,dur,[i,0],stim.vgrat(1),k0,stim.theta_g);
        end
end
if stim.disp == 1
    prompt = 'Press any number to start visualization of visual stimulus\n';
    start = input(prompt);
    figure
    if ~isempty(start)
        tmp = II{1};
        if isa(tmp,'plaid')
            tmp = generate_plaid(tmp);
        end
        for i=1:dur
            imagesc(squeeze(tmp(:,:,i)))
            drawnow
            pause(0.1)
        end
    end
end

%% motion-in-depth descriptors analysis - tuning curves
n_vel = length(param.prefVel);
n_orient = param.nOrient;
% [IIv,IIo] = size(II);
% ph_shift = 0;
% parameters{3,1} = ph_shift;
if (size(II,1)*size(II,2))>1
    e = zeros(4,1,1,n_orient,n_vel,size(param.phShift,2),size(II,1),size(II,2));
else
%    e=zeros(4,size(II{1},1),size(II{1},2),n_orient,n_vel,size(param.phShift,2),size(II,1),size(II,2));
end
TUN = 0;
% %ocular dominance is 0 or 1 'cause is monocular test
param.ocDom = 1; %force this value
for i=1:size(II,1)
    for j=1:size(II,2)
        
        disp([i j])
        %select input
        if isa(II{i,j},'plaid')
            tmp = generate_plaid(II{i,j});
            tmp = tmp(180:end-180,180:end-180,:);
            I = cat(4,zeros(size(tmp)),tmp);
        else
            I = cat(4,zeros(size(II{i,j})),II{i,j});
        end
        if ~debug
            [MT,EC21,EC22,EC1] = pop_flow_V1MT(I,param,TUN);
        else
            [MT,EC21,EC22,EC1,ES] = pop_flow_V1MT(I,param,TUN);
        end
        if (size(II,1)*size(II,2))>1
            if isa(II{i,j},'plaid')
                sze = size(MT);
                MT = squeeze(mean(mean(MT(sze(1)/2,sze(2)/2,:,:),1),2));
                EC21 = squeeze(mean(mean(EC21(sze(1)/2,sze(2)/2,:,:),1),2));
%                 EC22 = squeeze(mean(mean(EC22(sze(1)/2,sze(2)/2,:,:),1),2));
                EC1{1} = squeeze(mean(mean(EC1{1}(sze(1)/2,sze(2)/2,:,:),1),2));
                EC1{2} = squeeze(mean(mean(EC1{2}(sze(1)/2,sze(2)/2,:,:),1),2));
%                 EC1{3} = squeeze(mean(mean(EC1{3}(sze(1)/2,sze(2)/2,:,:),1),2));
%                 EC1{4} = squeeze(mean(mean(EC1{4}(sze(1)/2,sze(2)/2,:,:),1),2));
            else
                MT = squeeze(mean(mean(MT,1),2));
                EC21 = squeeze(mean(mean(EC21,1),2));
%                 EC22 = squeeze(mean(mean(EC22,1),2));
                EC1{1} = squeeze(mean(mean(EC1{1},1),2));
                EC1{2} = squeeze(mean(mean(EC1{2},1),2));
%                 EC1{3} = squeeze(mean(mean(EC1{3},1),2));
%                 EC1{4} = squeeze(mean(mean(EC1{4},1),2));
            end
        end

        if debug
            e(8,:,:,:,:,:,i,j) = ES(1,:,:,:,:,:);
            e(7,:,:,:,:,:,i,j) = ES(2,:,:,:,:,:);
            e(6,:,:,:,:,:,i,j) = ES(3,:,:,:,:,:);
            e(5,:,:,:,:,:,i,j) = ES(4,:,:,:,:,:);
        end
        e(4,:,:,:,:,:,i,j) = EC1{1};
        e(3,:,:,:,:,:,i,j) = EC1{2};
        e(2,:,:,:,:,:,i,j) = EC21; 
        e(1,:,:,:,:,:,i,j) = MT; 
        fprintf('Pop-activity, stimulus %d \n',size(II,2)*(i-1)+j);
        clear I
    end
end

e = squeeze(e);

varargout{1} = e;
varargout{2} = param;
%Save data in SIMULATIONS Directory
% path = 'SIMULATIONS/mono';
% OldFolder = cd;
% cd(path);
% save(file_name,'e','param','-v7.3')
% cd(OldFolder)
end

function argCheck(stimulus,fieldName)
j=1;
err = [];
for i=1:length(fieldName)
    if ~isfield(stimulus,fieldName{i})
        err(j) = i;
        j=j+1;
    end
end
if ~isempty(err)
    error('Following stimulus properties are not defined:\n %s %s',fieldName{err})
else
end
end
