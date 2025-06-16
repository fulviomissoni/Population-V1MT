function stim = initPlaidStimulus(truetheta,theta,vpld,contr,k,samples,k_gauss,varargin)
% arg.dur = dur;                              %aperture size      [pixs]
% arg.apert_rad = ceil(samples/2)+2;          %aperture size      [pixs]
% arg.truetheta = stim.truetheta(num_pld);    %true orientation   [rad]
% arg.vpld = stim.vel_stim(num_pld);          %velocity amplitude [pixs/frame]
% arg.k = [k0,k0];                            %spatial freq       [cycle/pix]
% arg.vgrat = stim.vgrat(num_pld,:);          %gratings vel       [pixs/frame]
% arg.theta_g = stim.theta_g(num_pld,:);      %gratings orient    [rad]
% arg.alpha = 0.5;                            %alpha channel for transparency
% arg.contrast                                %Contrast of two gratings (if
%                                             numeric then is used as
%                                             difference between the two
%                                             gratings
% arg.mode = stim.mode;                       %stimulus implementation algorithm
% arg.pl_type = pl_type;                      %plaid type
% arg.k_gauss = stim.k_gauss;                 %with k you can determine the size of the filter that will blur the image
%                                             %size = 1 / (k_gauss * k0)
% k:    spatial frequency of the gratings, if is scalar spatial frequency
% will be the same for both
%samples: size of the matrix in samples
%alpha (optional): transparency


if isempty(varargin)
    stim.alpha = 0.5;       %alpha channel for transparency
end

    if size(contr,2)<2
        c = 1/4+[-contr/4,contr/4];
    end
    if isscalar(k)
        k0 = [k,k];
    end
    stim.type = 'plaid';
    stim.truetheta =  truetheta(:); %true orientation
    stim.theta_g = theta; %true orientation
    stim.vel_stim = [vpld(:)];
    [x, y, z, c1] = ndgrid(stim.truetheta, stim.theta_g(:,1), stim.vel_stim, c(:,1));
    y = pagetranspose(y);   %pagetranspose serve per rendere le grid create con ndgrid nello stesso ordine dei meshgrid
    c1 = pagetranspose(c1);
    [stim.truetheta, y2, stim.vel_stim,c2] = ndgrid(stim.truetheta, stim.theta_g(:,2), stim.vel_stim, c(:,2));
    stim.truetheta = pagetranspose(stim.truetheta);
    y2 = pagetranspose(y2);
    stim.vel_stim = pagetranspose(stim.vel_stim);
    c2 = pagetranspose(c2);
    stim.truetheta = stim.truetheta(:); 
    stim.theta_g = stim.truetheta + [y(:) y2(:)];
    stim.contrast = [c1(:),c2(:)];
    stim.vel_stim = stim.vel_stim(:);
    if size(stim.theta_g,2)==2
        stim.ori = [cos( ...
            wrapTo2Pi(stim.theta_g(:,1)-stim.truetheta)),... 
            cos( ...
            wrapTo2Pi(stim.theta_g(:,2)-stim.truetheta))];
        stim.vgrat = [stim.ori(:,1).*stim.vel_stim, stim.ori(:,2).*stim.vel_stim];
        stim.vgrat = round(stim.vgrat,5,"decimals");
    else
        stim.vgrat = stim.vel_stim;
    end
    stim.dur = 43; %duration in frame
    stim.mode = 1;
    stim.disp = 0; %set to 1 to show visual stimulus in a figure
    stim.k = k0;
    stim.apert_rad = ceil(samples/2)+2;
    stim.vpld = vpld;
    stim.c = c; 
    stim.k_gauss = k_gauss;
end