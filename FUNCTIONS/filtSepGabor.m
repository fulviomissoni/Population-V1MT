function IF = filt_sep_gabor(I,fname,DC_thr)
% Compute the 2D filtering of the stereo input with 2D Gabor filters (FI_n)
% Note that 2D filtering is achieved as a combination of two 1D filtering
%           I       stereo image(input);
%           fname   name of .mat file that contains the 1D Gabor components
%                   needed to obtain 16 orientation channels (theta =
%                   0:pi/8:2*pi-pi/8);
%           DC_thr  threshold value. Default equal to "0";
%
% FUNDAMENTAL NOTE: the orientation value, in this function, are defined
% respect to the image reference system (centered at the top left). When
% you plot the shape filter the (that uses this frame of reference)
% orientation value is correct.

%Set DC threshold
if (nargin<3)
  DC_thr = 1e-10;
end


fltr = load(fname);
taps = fltr.taps;
F = fltr.F;
clear fltr;


n_orient = 16;

[nr,nc,n_frames] = size(I);
padsize = floor(taps/2); 
IF = zeros(nr,nc,n_orient,n_frames);
I = reshape(I,nr,nc,n_frames);

for ind = 1:n_frames


      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Horizontal and vertical %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%


      % 0

      F1Y = conv2b(I(:,:,ind),F(1,:).',3);
      even = conv2b(F1Y,F(2,:),3);
      odd = conv2b(F1Y,F(3,:),3);

      IF(:,:,1,ind) = complex(even,odd);
      
      %pi
      
      IF(:,:,9,ind) = complex(even,-odd);
      clear even odd F1Y;

      % pi/2

      F1X = conv2b(I(:,:,ind),F(1,:),3);
      even = conv2b(F1X,F(2,:).',3);
      odd = conv2b(F1X,F(3,:).',3);

      IF(:,:,5,ind) = complex(even,odd);

      % 11 pi/8
      IF(:,:,13,ind) = complex(even,-odd);

      clear even odd F1X;

%       % 11 pi/8
% 
%       F1X = conv2b(I(:,:,ind),F(1,:),3);
%       even = conv2b(F1X,F(2,:).',3);
%       odd = -conv2b(F1X,F(3,:).',3);


%       clear even odd F1X;


      %%%%%%%%%%%%
      % Diagonal %
      %%%%%%%%%%%%


      F4Y = conv2b(I(:,:,ind),F(4,:).',3);

      F4YF4X = conv2b(F4Y,F(4,:),3);
      F4YF5X = conv2b(F4Y,F(5,:),3);

      clear F4Y;

      F5Y = conv2b(I(:,:,ind),F(5,:).',3);

      F5YF4X = conv2b(F5Y,F(4,:),3);
      F5YF5X = conv2b(F5Y,F(5,:),3);

      clear F5Y;

      % pi/4

      even = F4YF4X - F5YF5X;
      odd = F5YF4X + F4YF5X;

      IF(:,:,3,ind) = complex(even,odd);

        % 10 pi/8

      even = F4YF4X - F5YF5X;
      odd = -(F5YF4X + F4YF5X);

      IF(:,:,11,ind) = complex(even,odd);

      % 3pi/4

      even = F4YF4X + F5YF5X;
      odd = F5YF4X - F4YF5X;

      IF(:,:,7,ind) = complex(even,odd);



        % 14pi/4

      even = F4YF4X + F5YF5X;
      odd = -(F5YF4X - F4YF5X);

      IF(:,:,15,ind) = complex(even,odd);

      clear even odd F4YF4X F4YF5X F5YF4X F5YF5X;


      %%%%%%%%%%%%%%%%%
      % 'In-betweens' %
      %%%%%%%%%%%%%%%%%


      F8Y = conv2b(I(:,:,ind),F(8,:).',3);

      F8YF6X = conv2b(F8Y,F(6,:),3);
      F8YF7X = conv2b(F8Y,F(7,:),3);

      clear F8Y;

      F9Y = conv2b(I(:,:,ind),F(9,:).',3);

      F9YF6X = conv2b(F9Y,F(6,:),3);
      F9YF7X = conv2b(F9Y,F(7,:),3);

      clear F9Y;

      % pi/8

      even = F8YF6X - F9YF7X;
      odd = F9YF6X + F8YF7X;

      IF(:,:,2,ind) = complex(even,odd);

      % 9pi/8

      even = F8YF6X - F9YF7X;
      odd = -(F9YF6X + F8YF7X);

      IF(:,:,10,ind) = complex(even,odd);

      % 7pi/8

      even = F8YF6X + F9YF7X;
      odd = F9YF6X - F8YF7X;

      IF(:,:,8,ind) = complex(even,odd);

      % 15pi/8

      even = F8YF6X + F9YF7X;
      odd = -(F9YF6X - F8YF7X);

      IF(:,:,16,ind) = complex(even,odd);

      clear even odd F8YF6X F8YF7X F9YF6X F9YF7X;


      F6Y = conv2b(I(:,:,ind),F(6,:).',3);

      F6YF8X = conv2b(F6Y,F(8,:),3);
      F6YF9X = conv2b(F6Y,F(9,:),3);

      clear F6Y;

      F7Y = conv2b(I(:,:,ind),F(7,:).',3);

      F7YF8X = conv2b(F7Y,F(8,:),3);
      F7YF9X = conv2b(F7Y,F(9,:),3);

      clear F7Y;

      % 3pi/8

      even = F6YF8X - F7YF9X;
      odd = F7YF8X + F6YF9X;

      IF(:,:,4,ind) = complex(even,odd);

      % 11pi/8

      even = F6YF8X - F7YF9X;
      odd = -(F7YF8X + F6YF9X);

      IF(:,:,12,ind) = complex(even,odd);

      % 5pi/8

      even = F6YF8X + F7YF9X;
      odd = F7YF8X - F6YF9X;

      IF(:,:,6,ind) = complex(even,odd);

      % 13pi/8

      even = F6YF8X + F7YF9X;
      odd = -(F7YF8X - F6YF9X);

      IF(:,:,14,ind) = complex(even,odd);

      clear even odd F6YF8X F6YF9X F7YF8X F7YF9X;

end

IF = reshape(IF,nr,nc,n_orient,n_frames);
invalid = (abs(real(IF))<DC_thr) | (abs(imag(IF))<DC_thr);
% IF(invalid) = NaN;
IF(invalid) = 0;

G_even = real(IF);
G_odd = imag(IF);
clear IF;
IF{1} = G_even;
IF{2} = G_odd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%