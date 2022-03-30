function p = plaid(arg)
% MOD = PLAID(THETA1, THETA2): creates a plaid object
% PLAID specifies geometry of plaid movements 
% NB: theta1 and theta2 are expressed in rad

% [p.C1,p.R1,p.C1S]=grating_geometry(grat.theta(1));
% [p.C2,p.R2,p.C2S]=grating_geometry(grat.theta(2));
if ~(length(arg.k)==2)
    error('Spatial frequency must be defined for both gratings!!')
end
if ~(length(arg.vgrat)==2)
    error('Velocity must be defined for both gratings!!')
end
if ~(length(arg.theta_g)==2)
    error('Velocity must be defined for both gratings!!')
end
p.apert_rad = arg.apert_rad;
p.dur = arg.dur;
p.truetheta = arg.truetheta;
p.vpld = arg.vpld;
p.k = arg.k;
p.vgrat = arg.vgrat;
p.theta_g = arg.theta_g;
p.alpha = arg.alpha;
p.c = arg.contrast;
p.pl_type = arg.pl_type;
p.mode = arg.mode;
% from structure to class...
p = class(p,'plaid');