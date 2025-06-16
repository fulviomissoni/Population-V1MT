function [varargout] = motionPopV1MT(param,stim)

samples = param.samples;
k0 = param.spat_freq;
%assuming same duration and type of stimuli
dur = stim(1).dur;
type = stim(1).type;
if ~isfield(stim,'type')
    error('Define type of the stimulus!! Options are: moving grats, plaids, or moving RDS')
end
field_names{1} = 'dur'; 
field_names{2} = 'vgrat'; 
% field_names{3} = 'theta_g';
% fieldNames{4} = 'truetheta'; fieldNames{5} = 'contrast_g';
% fieldNames{6} = 'mode'; fieldNames{7} = 'vel_stim';
argCheck(stim,field_names);
%% Stimulus definition 
stimuli = ["plaid","grat","RDS_moving","shift_grat","plaid_noise"];
if ~sum(matches(stimuli,type))
    error('Select from the possible stimuli:\n %s %s %s %s',stimuli{1},stimuli{2},stimuli{3},stimuli{4})
end
    
switch type
    case 'grat'
        II = sinGrating(samples,samples,dur,[0,0],stim.vgrat,k0,stim.theta_g);
    case 'plaid'
        % if ~isfield(stim,'pl_type')
        %     pl_type = 1;
        % else 
        %     pl_type = stim.pl_type;
        % end
        
        for num_pld = 1:length(stim)
            %plaid object
            % arg.dur = dur;                              %aperture size      [pixs]
            % arg.apert_rad = ceil(samples/2)+2;          %aperture size      [pixs]
            % arg.truetheta = stim.truetheta(num_pld);    %true orientation   [rad]
            % arg.vpld = stim.vel_stim(num_pld);          %velocity amplitude [pixs/frame]
            % arg.k = [k0,k0];                            %spatial freq       [cycle/pix]
            % arg.vgrat = stim.vgrat(num_pld,:);          %gratings vel       [pixs/frame]
            % arg.theta_g = stim.theta_g(num_pld,:);      %gratings orient    [rad]
            % % arg.alpha = 0.5;                            %alpha channel for transparency
            % % arg.contrast = stim.contrast_g(num_pld,:);  %Contrast of two gratings
            % arg.mode = stim.mode;                       %stimulus implementation algorithm
            % arg.pl_type = pl_type;                      %plaid type
            % arg.k_gauss = stim.k_gauss;                 %with k you can determine the size of the filter that will blur the image
            %                                             %size = 1 / (k_gauss * k0)
            %define plaid object
            II{num_pld} = plaid(stim(num_pld));
            %generate plaid stimulus
        end
    case 'RDS_moving'         
        for num_stim = 1:length(stim.truetheta)
                vx = stim.vgrat(num_stim,1).*cos(stim.truetheta(num_stim));
                vy = stim.vgrat(num_stim,2).*sin(stim.truetheta(num_stim));
                % scale_ind = 4; %do not remember what scale means, it should be related to different scales of the image from same size to bigger ones
                scale_ind = 1;
                II{num_stim} = movingRDS_MS(samples,samples,dur,scale_ind,vx, vy);
                % II{num_stim} = II{num_stim}(60:end-60,60:end-60,:);
        end
        % II = reshape(II,stim.stim_size);
    case 'shift_grat'
        II = sinGrating(samples,samples,dur,[0,0],stim.vgrat(1),k0,stim.theta_g);
        for i=1:floor(1/(2*k0))+1
           II(i+1) = sinGrating(samples,samples,dur,[i,0],stim.vgrat(1),k0,stim.theta_g);
        end

    case 'plaid_noise'
        if ~isfield(stim,'pl_type')
            pl_type = 1;
        else 
            pl_type = stim.pl_type;
        end
        
        for num_stim = 1:length(stim.truetheta)
            %plaid object
            arg.dur = dur;                              %aperture size      [pixs]
            arg.apert_rad = ceil(samples/2)+2;          %aperture size      [pixs]
            arg.truetheta = stim.truetheta(num_stim);    %true orientation   [rad]
            arg.vpld = stim.vel_stim(num_stim);          %velocity amplitude [pixs/frame]
            arg.k = [k0,k0];                            %spatial freq       [cycle/pix]
            arg.vgrat = stim.vgrat(num_stim,:);          %gratings vel       [pixs/frame]
            arg.theta_g = stim.theta_g(num_stim,:);      %gratings orient    [rad]
            arg.alpha = 0.5;                            %alpha channel for transparency
            arg.contrast = stim.contrast_g(num_stim,:);  %Contrast of two gratings
            arg.mode = stim.mode;                       %stimulus implementation algorithm
            arg.pl_type = pl_type;                      %plaid type
            arg.k_gauss = stim.k_gauss;                      %with k you can determine the size of the filter that will blur the image
                                                        %size = 1 / (k_gauss * k0)
            %define plaid object
            II_plaid{num_stim} = plaid(arg);
            %generate plaid stimulus
            I_plaid{num_stim} = generate_plaid(II_plaid{num_stim});

            %RDS_moving
            % scale_ind = 4; %do not remember what scale means, it should be related to different scales of the image from same size to bigger ones
            scale_ind = 1;
            vx_grat1 = stim.vgrat(num_stim,1).*cos(stim.theta_g(num_stim,1));
            vy_grat1 = stim.vgrat(num_stim,2).*sin(stim.theta_g(num_stim,1));
            I_RDS_grat1{num_stim} = movingRDS_MS(size(I_plaid{num_stim},1),size(I_plaid{num_stim},2), ...
                size(I_plaid{num_stim},3),scale_ind,vx_grat1, vy_grat1).*stim.sigma_noise;
            vx_grat2 = stim.vgrat(num_stim,1).*cos(stim.theta_g(num_stim,2));
            vy_grat2 = stim.vgrat(num_stim,2).*sin(stim.theta_g(num_stim,2));
            I_RDS_grat2{num_stim} = movingRDS_MS(size(I_plaid{num_stim},1),size(I_plaid{num_stim},2), ...
                size(I_plaid{num_stim},3),scale_ind,vx_grat2, vy_grat2).*stim.sigma_noise;
            % II{num_stim} = II{num_stim}(60:end-60,60:end-60,:);
            II{num_stim} = I_plaid{num_stim} + I_RDS_grat1{num_stim} + I_RDS_grat2{num_stim};
            
            II{num_stim} = (II{num_stim} - min(min(II{num_stim})))./ ...
                (max(max(II{num_stim})) - min(min(II{num_stim})));
        end

end
% if stim.disp == 1
%     prompt = 'Press any number to start visualization of visual stimulus\n';
%     start = input(prompt);
%     figure
%     if ~isempty(start)
%         tmp = II{1};
%         if isa(tmp,'plaid')
%             tmp = generate_plaid(tmp);
%         end
%         for i=1:dur
%             imagesc(squeeze(tmp(:,:,i)))
%             drawnow
%             pause(0.1)
%         end
%     end
% end

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
            I = II{i,j};
        end
        %population processing
        EC1 = popFlowV1MT(I,param);
        sze = size(EC1{1});
        if strcmp(type,"RDS_moving")
            %then we have to make th average to get to actual neuron
            %response
            e(1,:,:,i,j) = squeeze(mean(EC1{1}( ...
                ceil(param.filter_sample/2):param.filter_sample:end, ...
                ceil(param.filter_sample/2):param.filter_sample:end, ...
                :, ...
                :), ...
                [1,2], 'omitnan'));
            e(2,:,:,i,j) = squeeze(mean(EC1{2}( ...
                ceil(param.filter_sample/2):param.filter_sample:end, ...
                ceil(param.filter_sample/2):param.filter_sample:end, ...
                :, ...
                :), ...
                [1,2], 'omitnan'));
        else
            %in this case no average as the stimulus is without noise
            e(1,:,:,i,j) = squeeze(EC1{1}(ceil(sze(1)/2), ceil(sze(2)/2),:,:));
            e(2,:,:,i,j) = squeeze(EC1{2}(ceil(sze(1)/2), ceil(sze(2)/2),:,:));
        end

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