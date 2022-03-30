function surf_MID_tun(varargin)
%Surface plot of motion-in-depth of a sub-population of C21,C22 and MT
%cells. In the tuning analysis provided the number of orientation 
%channel, phase shift channel and v_pref is one (see parameters).
%It can use another population activity file in .mat extension identified 
%by its file name. In this file must be provided several information:
%   e:          multidimensional matrix that contains the population activity 
%               size: ncell x n_orient x n_ph x v_pref x vL_samples x vR_samples
%               where ncell is the number of cells, n_orient is the number 
%               of the orientation channels, n_ph is the number of the phase-shift channels, 
%               v_pref is the number of the velocity channels, and
%               vL_samples and vR_samples are the number of samples of the
%               velocity values of the moving sinusoidal gratings.
%   v:          velocity values vector
%   parameters: see pop_flowV1MT methods
%
%Fulvio Missoni

if nargin == 1
    file_name = varargin{1};
else
    file_name='MID_SurfTuning';
end
%load data
root='SIMULATIONS/MID-tuning/';
file_ext='.mat';
load([root file_name file_ext],'e','parameters','v');
%parameters
[n_cell,n_o,n_v,n_p,v1,v2]=size(e);
v_pref = parameters{4,1};
% psi = ["pi/4","0","-pi/4"];
psi = parameters{3,1};
% v1=v1(46:141);
cell = ['MT_';'C21';'C22';'C11';'C12';'C13';'C14'];
for i=1:n_cell
    for j=1:n_v
        for k=1:n_p
            my_txt = [cell(i,:),' v_c = ',num2str(v_pref(j)),' \Delta\psi = ',num2str(psi(k))];
            figure, imagesc(v,v,squeeze(e(i,1,j,k,:,:)))
%             colormap gray
            m = min(min(squeeze(e(i,1,j,k,:,:))));
            M = max(max(squeeze(e(i,1,j,k,:,:))));
            if i>2
                h=colorbar;
%                 set(h,'Ticks',[m,0,M],'TickLabels',[-1 0 1])
            end
            view(0,-90)
            title(my_txt)
            xlabel('v_L [px/frame]')
            ylabel('v_R [px/frame]')
        end
    end
end