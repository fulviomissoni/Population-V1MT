function mono_motion_simulation(param)

file_name='mono_motion_tun';
samples=param.samp;
if samples==11
    input_file = 'plaid43x43x43_vtrue0_2_k0_0_025';
end
if samples==43
    input_file = 'plaid43x43x43_vtrue_2_k0_0_063';
end
root='IMAGES/input/';
img_type='.mat';
load([root input_file img_type],'II');
%% motion-in-depth descriptors analysis - tuning curves
n_vel = length(param.prefVel);
n_orient = param.nOrient;
% ph_shift = 0;
% parameters{3,1} = ph_shift;
e=zeros(4,size(II{1},1),size(II{1},2),n_orient,n_vel,size(param.phShift,2),size(II,1)); 
TUN = 0;
%ocular dominance is 0 or 1 'cause is monocular test
param.ocDom = 0;
for i=1:size(II,1)
    %select input
    I = cat(4,II{i},zeros(size(II{i})));
    [MT,EC21,EC22,EC1] = pop_flow_V1MT(I,param,TUN);
    e(4,:,:,:,i) = squeeze(EC1{1});
    e(3,:,:,:,i) = squeeze(EC1{2});
    e(2,:,:,:,i) = squeeze(EC21); 
    e(1,:,:,:,i) = squeeze(MT); 
    fprintf('%d \n',i);
    clear I
end
%Save data in SIMULATIONS Directory
path = ['SIMULATIONS/mono'];
OldFolder = cd;
cd(path);
save(file_name,'e','param','-v7.3')
cd(OldFolder)
