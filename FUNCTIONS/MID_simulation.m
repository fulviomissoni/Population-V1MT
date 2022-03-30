function MID_simulation(param,debug)

if debug
    error('This method is under revision!!')
end
file_name='MID_SurfTuning';
ok=0;
if param.samples == 11
    input_file = 'sinGrating_11x11x11_k0_0_25';
    ok = 1;
end
if param.samples == 43
    input_file = 'sinGrating_43x43x43_k0_0_063';
    ok = 1;
end
if ok == 0
    error('Suitable stimulus file is not in folder!!')
end

root='IMAGES/input/';
img_type='.mat';
load([root input_file img_type],'II','v');
%% motion-in-depth descriptors analysis - tuning curves
n_orient = param.nOrient;
n_vel = length(param.prefVel);
% ph_shift = 0;
% parameters{3,1} = ph_shift;
e=zeros(7,n_orient,n_vel,size(param.phShift,2),size(II,1),size(II,1)); 
TUN = 1;

for ivL=1:size(II,1)
    tic
    for ivR=1:size(II,1)
        %select input
        I = cat(4,II{ivL},II{ivR});
        [MT,EC21,EC22,EC1] = pop_flow_V1MT(I,param,TUN);
        e(7,:,:,:,ivL,ivR) = squeeze(EC1{1});
        e(6,:,:,:,ivL,ivR) = squeeze(EC1{2});
        e(5,:,:,:,ivL,ivR) = squeeze(EC1{3});
        e(4,:,:,:,ivL,ivR) = squeeze(EC1{4});
        e(3,:,:,:,ivL,ivR) = squeeze(EC21);
        e(2,:,:,:,ivL,ivR) = squeeze(EC22); 
        e(1,:,:,:,ivL,ivR) = squeeze(MT); 
        fprintf('%d %d\n',ivL,ivR);
        clear I
    end
    toc
end
%Save data in SIMULATIONS Directory
path = ['SIMULATIONS/MID-tuning'];
OldFolder = cd;
cd(path);
save(file_name,'e','param','v')
cd(OldFolder)
