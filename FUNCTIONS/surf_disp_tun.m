function surf_disp_tun(varargin)
%Surface plot of disparity tuning of a sub-population of disparity
%detectors (C21, C22). In the tuning analysis provided the number of orientation 
%channels are eight (span from 0 to pi) and the number of phase-shift values 
%are seven.
%It can use another population activity file in .mat extensione identified 
%by its file name. In this file must be provided several information:
%   e:          multidimensional matrix that contains the population activity 
%               size: ncell x n_orient x n_ph x d_samples x d_samples x
%               Nseed where ncell is the number of cells, n_orient is the
%               number of the orientation channels, phase_num is the number of
%               the phase-shift channels, d_samples is the number of samples of
%               the disparity vector and Nseed is the number of seeds of the
%               pseudo-random generator.
%   d:          disparity values vector
%   c_or:       string that contains the labels of each orientation channel
%   c_ph:  string that contains the labels of each phase-shift channel

if nargin == 1
    file_name = varargin{1};
else
    file_name='RDS_SurfTuning2';
end
%load data
root='SIMULATIONS\disparity-tuning\';
file_ext='.mat';
load([root file_name file_ext],'e','d','c_ph','c_or');
%other parameters
[ncell, n_orient, n_ph, d_samples, d_samples, Nseed]=size(e);
%string for plot
for c=1:ncell
    h=figure;
    for o=1:n_orient
        for ph=1:n_ph
            h1=subplot(n_orient,n_ph,ph+n_ph*(o-1));
            sumC = zeros(length(d),length(d));
            etmp = squeeze(e(c,o,ph,:,:,:));
            for i=1:Nseed
                sumC = sumC + squeeze(etmp(:,:,i));
            end
            sumC=sumC./Nseed;
            imagesc(d,d,sumC)
            set(gca,'PlotBoxAspectRatio',[1 0.45 1])
            hold on
            [m,indY,indX] = max2D(sumC);
            plot(xlim,[0 0],'k')
            plot([0 0],ylim,'k')
            % Mark empirical maximum of actual 2D disparity tuning surface:
            plot(d(indX),d(indY),'wx','linew',2,'markersize',10)
            view(0,-90)
            if o==1
                phasestr = c_ph{ph};
                text(0.5,1.42,['\Delta\phi = ' phasestr ],'FontSize',11,'units','norm','horizontalal','cent','verti','middle')
                if ph==1
                    text(-1,1.4,'\phi_0=\pi/2','units','norm','HorizontalAlignment','left','VerticalAlignment','middle','edgecolor','k')
                end
            end
            if ph==1
                orientstr = c_or{o};
                text(-0.58,0.5,['\theta = ' orientstr ],'FontSize',11,'rot',90,'units','norm','horizontalal','cent','verti','middle')
            end
            colormap jet
            axis off

        end
    end
    set(h,'position',[300 500 1300 850])
    h=colorbar;
    set(h,'pos',[0.93 0.1 0.015 0.85],'fonts',11)
end
end
%%%%%   %%%%%   %%%%%   %%%%%   %%%%%   %%%%%   %%%%%   %%%%%   %%%%%   %%%
function [m,indY,indX] = max2D(A)
    if length(size(A))>2||length(size(A))<2
        error('size of the matrix must be equal to 2!!')
    end
    [mtmp,indY] = max(A);
    [m,indX] = max(mtmp);
    indY = indY(indX);
end