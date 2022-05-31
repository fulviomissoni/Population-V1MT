function [varargout] = motionPopV1MT(param,stim)

samples = param.samples;
k0 = param.spat_freq;
dur = stim.dur;
if ~isfield(stim,'type')
    error('Define type of the stimulus!! Options are: moving grats, plaids, or moving RDS')
end
field_names{1} = 'dur'; 
field_names{2} = 'vgrat'; field_names{3} = 'theta_g';
% fieldNames{4} = 'truetheta'; fieldNames{5} = 'contrast_g';
% fieldNames{6} = 'mode'; fieldNames{7} = 'vel_stim';
argCheck(stim,field_names);
%% Stimulus definition 
stimuli = ["plaid","grat","RDS_moving","shift_grat"];
if ~sum(matches(stimuli,stim.type))
    error('Select from the possible stimuli:\n %s %s %s %s',stimuli{1},stimuli{2},stimuli{3},stimuli{4})
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
                scale_ind = 4;
                II{num_stim} = movingRDS_MS(samples,samples,dur,scale_ind,vx, vy);
                II{num_stim} = II{num_stim}(60:end-60,60:end-60,:);
        end
        II = reshape(II,stim.stim_size);
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
n_vel = length(param.pref_vel);
n_orient = param.n_orient;

%allocate memory for pop_response of only one cell (centered in image
%center)
e = zeros(2,n_orient,n_vel,size(II,1),size(II,2));

% %ocular dominance is 0 or 1 'cause is monocular test
for i=1:size(II,1)
    for j=1:size(II,2)
        
        disp([i j])
        %select input
        if isa(II{i,j},'plaid')
            tmp = generate_plaid(II{i,j});
%             tmp = tmp(180:end-180,180:end-180,:);
            I = tmp;
        else
            I = zeros(size(II{i,j}));
        end
        %population processing
        EC1 = popFlowV1MT(I,param);
        sze = size(EC1{1});

        e(1,:,:,i,j) = squeeze(EC1{1}(ceil(sze(1)/2), ceil(sze(2)/2),:,:));
        e(2,:,:,i,j) = squeeze(EC1{2}(ceil(sze(1)/2), ceil(sze(2)/2),:,:));
        fprintf('Pop-activity, stimulus %d \n',size(II,2)*(i-1)+j);
        clear I
    end
end

e = squeeze(e);

varargout{1} = e;
varargout{2} = param;

end
%% LOCAL FUNCTIONS
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
    end
end