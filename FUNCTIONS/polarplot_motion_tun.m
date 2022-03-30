function polarplot_motion_tun(varargin)
%
%Fulvio Missoni (2021)

if nargin == 1
    file_name = varargin{1};
else
    file_name='mono_motion_tun';
end

root='SIMULATIONS/mono/';
file_ext='.mat';
load([root file_name file_ext],'e','param');
[n_cell,sy,sx,n_orient,n_ph,N_vel,N_dir] = size(e);
thetatrue = 0:pi/8:(2*pi-pi/8);
mytitle=["MT";"C21";"C11";"C12"];
A = zeros(N_dir,1);
for j=1:n_orient
for i=1:1
    tmp = squeeze(e(i,21,21,j,1,1,:));
    figure, polarplot(thetatrue,tmp)
    title(mytitle(i))
end
end