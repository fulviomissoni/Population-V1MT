function II = generate_plaid(pld,varargin)
    if ~size(varargin,1)
        disp = 0;
    else
        if numel(varargin)>1
            if isempty(varargin{1})
                disp=0;
            else
                disp = varargin{1};
            end
        end
    end

    for j=1:size(pld.vgrat,1)
        II = simulate(pld);
    %         set(gcf,'Name','test','pos',[1 620 200 200])
        if disp
            figure
            set(gcf,'color','k')
            set(gca,'vis','off')
            for i=1:pld.dur 
                imagesc(II(:,:,i))
                view(0,-90)
                pause(0.05)
            end
        end
    end
end

function varargout = simulate(pl)

switch(pl.mode)
    case 1
        varargout{1} = analytic_express(pl);
    case 2
        [varargout{1}, varargout{2},varargout{3}]=im_rotate_mode(pl);
end
end
        
function II = analytic_express(pl)
    
% From plaid to grating movements
%     vtrue = pl.vpld*[cos(pl.truetheta); sin(pl.truetheta)]; %pixels/frame  this is true plaid velocity
%     m1 = C1*vtrue;   % degs/s  this is grating 1 velocity
%     m1n = sqrt(sum(m1.^2));
% 
%     m2 = C2*vtrue;   % degs/s  this is grating 2 velocity
%     m2n = sqrt(sum(m2.^2));
% 
%     omega1 = 2*pi*m1n; 
%     omega2 = 2*pi*m2n;
    [xx,yy]=meshgrid([-pl.apert_rad:pl.apert_rad],...
                     [-pl.apert_rad:pl.apert_rad]);
    ap_mask = xx.^2+yy.^2>pl.apert_rad.^2;

    [C1,R1,C1S]=grating_geometry(pl.theta_g(1));
    [C2,R2,C2S]=grating_geometry(pl.theta_g(2));

    xxr1 = R1*[xx(:)';yy(:)'];
    xxr2 = R2*[xx(:)';yy(:)'];

% simulates plaid movements
    c1 = pl.c(1);
    c2 = pl.c(2);
    i=1;
    for t=0:pl.dur
        pR1 = c1*(1+sign(cos(2*pi*pl.k(1)*xxr1(1,:)-(2*pi*pl.vgrat(1)*pl.k(1))*t)));

        pR2 = c2*(1+sign(cos(2*pi*pl.k(2)*xxr2(1,:)-(2*pi*pl.vgrat(2)*pl.k(2))*t)));

        pRe = reshape((pR1.*pl.alpha+(1-pl.alpha).*pR2),length(xx),length(yy));
        pRe = imfilter(pRe,fspecial('gaussian',ceil(length(xx)/3)));
        pRe(ap_mask)=ceil(mean([c1,c2]));
        i=i+1;
        II(:,:,i) = pRe;
    end
end
function [pRe,g1_rot,g2_rot]= im_rotate_mode(pl)
    warning('This modality of plaid generation has inexact velocity of plaid and gratings. This is due to the algorithm that is not suited for under-pixel movements.')
    prompt = 'Would you continue? [Y/n]\n';
    str = input(prompt,'s');
    if strcmp(str,'n')||strcmp(str,'N')
        error('The execution is stopped!! Set "mode" property of stim variable to 1')
    end
%     omega1 = 2*pi*m1n; 
%     omega2 = 2*pi*m2n;

%     nCycle =round(dim_x/pixPerCycle);
    nCycle = round(2*pl.apert_rad*pl.k(1));
    %recompute dim_x
    line_width(1) = floor(2*pl.apert_rad/(2*nCycle));
    Lmin = 0.5-pl.c(1); Lmax=0.5+pl.c(1);
    dim(1) = nCycle*line_width(1)*2;
    nCycle(2) = round(2*pl.apert_rad*pl.k(2));
    line_width(2) = floor(2*pl.apert_rad/(2*nCycle(2)));
    dim(2) = nCycle(2)*line_width(2)*2;
    dim = min(dim);
    g1 = repmat([Lmin*ones(dim,line_width(1)) Lmax*ones(dim,line_width(1))],1, nCycle(1));
    g1pad = padarray(g1,line_width*2,'circular','both');
    tmp = imrotate(g1pad,rad2deg(pl.theta_g(1)),'bicubic','loose');
    indStart = floor(size(tmp,1)/2-dim/2); indEnd = ceil(size(tmp,1)/2+dim/2);
    g1_rot = zeros(dim, dim, pl.dur);
    tmp(tmp<0) = 0;         tmp(tmp>1) = 1;
    g1_rot(:,:,1) = tmp(indStart:indEnd-1,indStart:indEnd-1);
%     g1 = g1(1:dim,1:dim);

    %recompute dim_x
%     dim_x = nCycle*pixPerCycle;
    Lmin = 0.5-pl.c(2); Lmax=0.5+pl.c(2);    
%     g2 = repmat([Lmin*ones(dim,line_width(2)) Lmax*ones(dim,line_width(2))],1, nCycle(2));
    g2 = repmat([Lmin*ones(dim,line_width(2)) Lmax*ones(dim,line_width(2))],1, nCycle(2));
    g2pad = padarray(g2,line_width*2,'circular','both');
%     g2 = g2(1:dim,1:dim);
    g2_rot = zeros(dim, dim, pl.dur);
    tmp = imrotate(g2pad,rad2deg(pl.theta_g(2)),'bicubic','loose');
    indStart = floor(size(tmp,1)/2-dim/2); indEnd = ceil(size(tmp,1)/2+dim/2);
    tmp(tmp<0) = 0;         tmp(tmp>1) = 1;
    g2_rot(:,:,1) = tmp(indStart:indEnd-1,indStart:indEnd-1);
    pRe = g1_rot;
%     g1_rot(:,:,1) = imrotate(g1,rad2deg(pl.theta_g(1)),'bicubic','crop');
    
%     g2_rot(:,:,1) = imrotate(g2,rad2deg(pl.theta_g(2)),'bicubic','crop');
    pRe(:,:,1) = reshape((g1_rot(:,:,1).*pl.alpha+(1-pl.alpha).*g2_rot(:,:,1)),dim,dim);
    [xx,yy] = meshgrid(1:dim,1:dim);
    ap_mask = (xx-ceil(dim/2)).^2 + (yy-ceil(dim/2)).^2 > ceil(dim/2)^2;
    ap_mask = repmat(ap_mask,[1 1 pl.dur]);
    %static gratings
    if ~pl.vgrat(1)
        return;
    end
    for j = 2:pl.dur           
        if ( floor(pl.vgrat(1)) == pl.vgrat(1) ) && ( floor(pl.vgrat(2)) == pl.vgrat(2) )
            g1 = circshift(g1, pl.vgrat(1), 2);
            g2 = circshift(g2, pl.vgrat(2), 2);
         else
            g1tmp = imtranslate(g1pad,[pl.vgrat(1),0],'bicubic','OutputView','full');
            g2tmp = imtranslate(g2pad,[pl.vgrat(2),0],'bicubic','OutputView','full');   
%             w = (size(g1tmp,2)-size(g1,2));
%             w(2) = (size(g2tmp,2)-size(g2,2));
            indStart = floor(size(g1tmp,1)/2-dim/2); indEnd = ceil(size(g1tmp,1)/2+dim/2);
            tmp = g1tmp(indStart:indEnd-1,indStart:indEnd-1); 
%             g1tmp(:,1:w(1)) = g1tmp(:,size(g1,2)+1:end);
%             g2tmp(:,1:w(2)) = g2tmp(:,size(g2,2)+1:end);
%             g1tmp(:,size(g1,2)+1:end)=[];
%             g2tmp(:,size(g2,2)+1:end)=[];
%             g1tmp(g1tmp>0.75)=1; g1tmp(g1tmp<=0.4)=0;
%             g2tmp(g2tmp>0.75)=1; g2tmp(g2tmp<=0.4)=0;
            g1 = tmp;  
            g1(g1<0) = 0;       g1(g1>1) = 1;
            indStart = floor(size(g2tmp,1)/2-dim/2); indEnd = ceil(size(g2tmp,1)/2+dim/2);
            tmp = g2tmp(indStart:indEnd-1,indStart:indEnd-1); 
            g2 = tmp;
            g2(g2<0) = 0;       g2(g2>1) = 1;
%             g1(g1<=0) = 0; g2(g2<=0) = 0;
%             g1 = g1tmp; g2 = g2tmp;
        end
        g1pad = padarray(g1,line_width*2,'circular','both');
        g2pad = padarray(g2,line_width*2,'circular','both');
        tmp = imrotate(g1pad,rad2deg(pl.theta_g(1)),'bicubic','crop');
        indStart = floor(size(tmp,1)/2-dim/2); indEnd = ceil(size(tmp,1)/2+dim/2);
        tmp(tmp<0) = 0;         tmp(tmp>1) = 1;
        g1_rot(:,:,j) = tmp(indStart:indEnd-1,indStart:indEnd-1); 
        
        tmp = imrotate(g2pad,rad2deg(pl.theta_g(2)),'bicubic','crop');
        indStart = floor(size(tmp,1)/2-dim/2); indEnd = ceil(size(tmp,1)/2+dim/2);
        tmp(tmp<0) = 0;         tmp(tmp>1) = 1;
        g2_rot(:,:,j) = tmp(indStart:indEnd-1,indStart:indEnd-1);
        pRe(:,:,j) = reshape((g1_rot(:,:,j).*pl.alpha+(1-pl.alpha).*g2_rot(:,:,j)),dim,dim);
    end
    pRe(ap_mask) = 0.5;
%     if length(speed) > 1
%         g1_rot(circle_index) = 0;
%         g2_rot(circle_index) = 0;
%     end
end

function [C,R,CS]=grating_geometry(theta)
% Define geometry of plaid movements
% theta = -theta;
% m = C*v
% theta =  ( pi/2 - theta ) ;
C = [cos(theta)^2 sin(theta)*cos(theta); ...
      sin(theta)*cos(theta) sin(theta)^2];


R = [cos(theta) -sin(theta);...
      sin(theta) cos(theta)]';
  
  
% mS = CS*v
CS = [cos(theta) sin(theta)];
end