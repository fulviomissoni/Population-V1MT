function plot_disp_tun(varargin)
%Horizontal disparity tuning curves plot of a sub-population of complex- and
%MT-like neurons varying the phase-shift channel.
%The orientation channel is fixed by the method (\theta = 0°) while the
%velocity channel is fixed by the user (only one value!).
%The function plot only the data stored in SIMULATIONS folder. By default the
%'disp_tuning_curve_norm.mat' file is plotted but it is possible to specify
%another file name(varargin{1}). This file must be contain a
%multidimensional matrix of name "e" and of size n_cell x n_orient x n_ph x N_d x N_seed
%where n_cell is the number of cells, n_orient is the number of the
%orientation channel (equal to 1), n_ph is the number of the phase shift
%channels, N_d is the number of samples of the disparity values vector and
%N_seed is the number of seeds of the pseudo-random generator
%
%Fulvio Missoni (2021)

if nargin == 1
    file_name = varargin{1};
else
    file_name='disp_tuning_curve_norm';
end

root='SIMULATIONS/disparity-tuning/';
file_ext='.mat';
load([root file_name file_ext],'e','d','parameters');
[n_cell,n_orient,n_ph,N_d,N_seed] = size(e);
%Control
if ~(n_orient == 1)
    error('The number of the orientation channels is fixed (only one)!')
end

if n_orient == 1
    mytitle=["MT";"C21";"C22";"C11";"C12";"C13";"C14"];
    for i=1:n_cell
        tmp = squeeze(e(i,:,:,:,:));
        tmp = mean(squeeze(tmp),length(size(tmp)));
        figure, plot(d,tmp')
%         ylim([0,3])
        title(mytitle(i))
    end
end