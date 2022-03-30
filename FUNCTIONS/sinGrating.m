function [II, n_frames, k0,t]=sinGrating(sx,sy,N,shift,v,k0,theta)
%   function to generate sinusoidal grating signal moving in temporal
%   domain -> each element of II contain sinuoidal grating at different
%   velocity
%
%   sx:     x dimension of input
%   sy:     y dimension of input
%   N:      number of taps
%   shift:  1x2 shift vector [dy,dx]
%   v:      velocity of the stimulus
%   theta:  orientation of the motion of the stimulus 
%           (default: pi/2)
%   k0:     frequency of the stimulus (default k0=.063 cycle/degree)

%   T (sampling step) -> fixed
%   
%   FULVIO MISSONI

if nargin == 5
    k0=.063;
end
if (nargin == 5)||(nargin == 6)
    theta=pi/2;
end
vt = length(v);
n_frames=N;
th = length(theta);
t=0:n_frames-1; %with this implementation each unit of the vector is a sample
II = cell(length(v),th); 
% II = cell(vt*th,1);
for f=1:length(v)
    for g=1:th
        II{f,g} = zeros(sy,sx,n_frames);
    end
end


% %Call  C++ routine
% II = MyMEXCode(length(v),v,length(theta),theta,k0,sx,sy,n_frames,t); ->
% the MEX function is not faster in all the cases (in general, matlab is
% faster)

%This is the fastest implementation
[y,x,t] = meshgrid(1:sy,1:sx,t);

% theta_stim = (theta-pi/2)-pi/2; %band orientation
% theta_stim = (theta-pi/2)+pi; %vector orientation
theta_stim = theta+pi/2;
% [v,theta_stim] = meshgrid(v,theta_stim);
% v = v(:); theta_stim = theta_stim(:);
noisetmp = (rand(sy,sx)*2-1);
for i=1:N
   noise(:,:,i) =  circshift(noisetmp,2,1);
end
for f=1:(vt) %for each velocities and orientation
    for g=1:th %for each orientation
        II{f,g} = 1 + sign(cos(2*pi*k0*((x+shift(2))*cos(theta_stim(g)) + (y+shift(1))*sin(theta_stim(g))) - (2*pi*k0)*v(f)*t)) + noise*0;
%         II{f,g} = 1+ noise*0.3;
    end
end
end
