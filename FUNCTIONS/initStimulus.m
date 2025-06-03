function stim = initStimulus(truetheta,theta,vpld,varargin)
    if numel(varargin)==2
        c(1) = varargin{1};
        c(2) = varargin{2};
    end
    if numel(varargin)==1
        c = 0.25+[-varargin{1}/4,varargin{1}/4];
    end
    stim.type = 'plaid';
    stim.truetheta =  truetheta(:); %true orientation
    stim.theta_g = [theta]; %true orientation
    stim.vel_stim = [vpld(:)];
    [x, y, z, c1] = ndgrid(stim.truetheta, stim.theta_g(:,1), stim.vel_stim, c(1));
    y = pagetranspose(y);   %pagetranspose serve per rendere le grid create con ndgrid nello stesso ordine dei meshgrid
    c1 = pagetranspose(c1);
    [stim.truetheta, y2, stim.vel_stim,c2] = ndgrid(stim.truetheta, stim.theta_g(:,2), stim.vel_stim, c(2));
    stim.truetheta = pagetranspose(stim.truetheta);
    y2 = pagetranspose(y2);
    stim.vel_stim = pagetranspose(stim.vel_stim);
    c2 = pagetranspose(c2);
    stim.truetheta = stim.truetheta(:); 
    stim.theta_g = stim.truetheta + [y(:) y2(:)];
    stim.contrast_g = [c1(:),c2(:)];
    stim.vel_stim = stim.vel_stim(:);
    if size(stim.theta_g,2)==2
        stim.ori = [cos(stim.theta_g(:,1)-stim.truetheta), cos(stim.theta_g(:,2)-stim.truetheta)];
        stim.vgrat = [stim.ori(:,1).*stim.vel_stim, stim.ori(:,2).*stim.vel_stim];
        stim.vgrat = round(stim.vgrat,5,"decimals");
    else
        stim.vgrat = stim.vel_stim;
    end
    stim.dur = 43; %duration in frame
    stim.mode = 1;
    stim.disp = 0; %set to 1 to show visual stimulus in a figure
end