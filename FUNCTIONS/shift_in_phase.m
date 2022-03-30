function [G_even, G_odd] = shift_in_phase(F_real,F_imag,ph_shift)
%Computes of response of neurons with different phase-shift values
%
%   F_real, F_imag      are, respectively, the real part and imaginary 
%                       part of the response of different orientation channels.
%
%   ph_shift            No x Np matrix where No is the number of the
%                       orientation channels and Np is the number of ph_shift 
%                       channels. phase-shift channels. Each row contains 
%                       the phase-shift values repeated for each orientation channel.
%
%

%% size of channels and size of the receptive fields
[nr,nc,n_frame,or,vv] = size(F_real);

% ALLOCATE MEMORY
G_even_tmp = reshape(F_real,nr*nc,n_frame*or*vv);
G_odd_tmp = reshape(F_imag,nr*nc,n_frame*or*vv);
ph_ind = 1;
% COMPUTE PHASE SHIFT
for ph=ph_shift(1,:)
    G_even(:,:,ph_ind) = G_even_tmp*cos(ph)-G_odd_tmp*sin(ph);
    G_odd(:,:,ph_ind) = G_odd_tmp*cos(ph)+G_even_tmp*sin(ph);
    ph_ind = ph_ind+1;
end
% RESHAPE OUTPUT RESULT
G_even = reshape(G_even,nr,nc,n_frame,or,vv,ph_ind-1);
G_odd = reshape(G_odd,nr,nc,n_frame,or,vv,ph_ind-1);
end
