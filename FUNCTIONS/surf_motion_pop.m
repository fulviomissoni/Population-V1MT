function surf_motion_pop(e,param)

if numel(size(e)) == 2
    [n_or,n_vel] = size(e);
    e = reshape(e,1,n_or,n_vel);
end
e = permute(e,[2 3 1]);
[n_or,n_vel,neuron]= size(e);

mytitle = cell(neuron);
mytext = cell(8);
mytext{1} = "MT";
mytext{2} = "C22";
mytext{3} = "C14";
mytext{4} = "C13";
mytext{5} = "S8";
mytext{6} = "S7";
mytext{7} = "S6";
mytext{8} = "S5";
for i = 1:neuron
    mytitle{i} = mytext{i};
end

 
if n_or==8
    theta_pref = linspace(0,pi,n_or+1);
    e = cat(1,e,flip(e(1,:,:),3));
end
if n_or==16
    theta_pref = linspace(0,2*pi,n_or+1);
    e = cat(1,e,e(1,:,:));
end
[s,t] = meshgrid(param.prefVel,theta_pref);
[xq,yq] = pol2cart(t,s);

for j=1:neuron
    tmp = squeeze(e(:,:,j));
    figure,
    surf(xq,yq,tmp,'EdgeColor','none','FaceColor','interp')
    axis off
    view(0,90)
    colorbar
    title(mytitle{j})
end
