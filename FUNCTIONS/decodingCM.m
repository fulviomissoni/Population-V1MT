function [vx,vy] = decodingCM(pop_resp_V1MT,param)

theta_cell_OUT = 0:pi/param.n_orient:pi-pi/param.n_orient;
[xx,tt] = meshgrid(param.pref_vel,theta_cell_OUT);
tt = tt + pi*(xx<0);
xx = abs(xx);
subPopResp = pop_resp_V1MT;
M = sum(sum(subPopResp));
%centre of mass decoding
vx = squeeze(sum(sum(subPopResp.*(xx.*cos(tt)),1),2)./M);
vy = squeeze(sum(sum(subPopResp.*(xx.*sin(tt)),1),2)./M);

end