function IF = filt_gabor_space2D(II,fname,flag,DC_thr)
% Compute the 2D filtering of the stereo input with 2D Gabor filters (FI_n)
% Note that 2D filtering is achieved as a combination of two 1D filtering
%           I       stereo image(input);
%           fname   name of .mat file that contains the 1D Gabor components
%                   needed to obtain 8 orientation channels (theta =
%                   0:pi/8:pi-pi/8);
%           flag    is set to "1" then is considered only the vertical
%                   orientation (theta = 0 rad);
%           DC_thr  threshold value. Default equal to "0";
%
% FUNDAMENTAL NOTE: the orientation value, in this function, are defined
% respect to the image reference system (centered at the top left). When
% you plot the shape filter the (that uses this frame of reference)
% orientation value is correct.

if (nargin<4)
    DC_thr =0;
end

%LOAD FILTER
fltr = load(fname); 

F=fltr.F;
clear fltr

%MEMORY ALLOCATION
[nr, nc, n_frames, side] = size(II);
taps = size(F,2); %numbers of taps of the filter

if ~flag
    n_orient = 8;
    IF{1} = zeros(nr,nc,n_orient,n_frames,side);
    IF{2} = zeros(nr,nc,n_orient,n_frames,side);
end
if flag
    IF = zeros(nr,nc,1,n_frames,side);
end

%CONVOLUTION COMPUTING
%for each sides
for s_ind=1:side
    %for each frames
    for frame=1:n_frames
        %compute the 2D response
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Horizontal and vertical %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%

        I=II(:,:,frame,s_ind);
        % 0 rad
        F1Y = conv2b(I,F(1,:).',3);
        even = conv2b(F1Y,F(2,:),3);
        odd = conv2b(F1Y,F(3,:),3);
        IF{1}(:,:,1,frame,s_ind) = even;
        IF{2}(:,:,1,frame,s_ind) = odd;
        clear even odd F1Y;
        
        if ~flag
            % pi/2 rad
            F1X = conv2b(I,F(1,:),3);
            even = conv2b(F1X,F(2,:).',3);
            odd = conv2b(F1X,F(3,:).',3);
            IF{1}(:,:,5,frame,s_ind) = even;
            IF{2}(:,:,5,frame,s_ind) = odd;
            clear even odd F1X;

            %%%%%%%%%%%%
            % Diagonal %
            %%%%%%%%%%%%

            F4Y = conv2b(I,F(4,:).',3);
            F4YF4X = conv2b(F4Y,F(4,:),3);
            F4YF5X = conv2b(F4Y,F(5,:),3);
            clear F4Y

            F5Y = conv2b(I,F(5,:).',3);
            F5YF4X = conv2b(F5Y,F(4,:),3);
            F5YF5X = conv2b(F5Y,F(5,:),3);
            clear F5Y;

            % pi/4
            even = F4YF4X - F5YF5X;
            odd = F5YF4X + F4YF5X;
            IF{1}(:,:,3,frame,s_ind) = even;
            IF{2}(:,:,3,frame,s_ind) = odd;

            % 3pi/4
             even = F4YF4X + F5YF5X;
             odd = F5YF4X - F4YF5X;
             IF{1}(:,:,7,frame,s_ind) = even;
             IF{2}(:,:,7,frame,s_ind) = odd;

            clear even odd F4YF4X F4YF5X F5YF4X F5YF5X;
            %%%%%%%%%%%%%%%%%
            % 'In-betweens' %
            %%%%%%%%%%%%%%%%%

            % pi/8
            F8Y = conv2b(I,F(8,:).',3);

            F8YF6X = conv2b(F8Y,F(6,:),3);
            F8YF7X = conv2b(F8Y,F(7,:),3);

            clear F8Y;

            F9Y = conv2b(I,F(9,:).',3);
            F9YF6X = conv2b(F9Y,F(6,:),3);
            F9YF7X = conv2b(F9Y,F(7,:),3);

            clear F9Y;

            % pi/8
            even = F8YF6X - F9YF7X;
            odd = F9YF6X + F8YF7X;

            IF{1}(:,:,2,frame,s_ind) = even;
            IF{2}(:,:,2,frame,s_ind) = odd;

            % 7pi/8
            even = F8YF6X + F9YF7X;
            odd = F9YF6X - F8YF7X;

            IF{1}(:,:,8,frame,s_ind) = even;
            IF{2}(:,:,8,frame,s_ind) = odd;

            clear even odd F8YF6X F8YF7X F9YF6X F9YF7X;        

            % 3pi/8
            F6Y = conv2b(I,F(6,:).',3);

            F6YF8X = conv2b(F6Y,F(8,:),3);
            F6YF9X = conv2b(F6Y,F(9,:),3);

            clear F6Y;

            F7Y = conv2b(I,F(7,:).',3);

            F7YF8X = conv2b(F7Y,F(8,:),3);
            F7YF9X = conv2b(F7Y,F(9,:),3);
            clear F7Y;

            % 3pi/8
            even = F6YF8X - F7YF9X;
            odd = F7YF8X + F6YF9X;

            IF{1}(:,:,4,frame,s_ind) = even;
            IF{2}(:,:,4,frame,s_ind) = odd;

            % 5pi/8
            even = F6YF8X + F7YF9X;
            odd = F7YF8X - F6YF9X;

            IF{1}(:,:,6,frame,s_ind) = even;
            IF{2}(:,:,6,frame,s_ind) = odd;

            clear even odd F6YF8X F6YF9X F7YF8X F7YF9X;
        end
    end
end

%THRESHOLDING
invalid = (abs(IF{1})<DC_thr) | (abs(IF{2})<DC_thr);
IF{1}(invalid) = NaN;
IF{2}(invalid) = NaN;

