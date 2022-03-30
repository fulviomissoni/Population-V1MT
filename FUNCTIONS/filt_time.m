function F_new = filt_time(F,conv_type,filter_type,v,k0)
% convolution with temporal component of the spatio-temporal filter:
%       F           convolution component with spatial filters
%       conv_type   type of convolution ('full','same','valid') -> see conv2
%                   help
%       filter_type temporal filter choice ('gabor','exp_decay','adelson-bergen') 
%                   (see 'technical report', Fulvio Missoni)
%       v           preferred velocity (NOTE: of the future detectors)
%       k0          gabor spatial frequency in cycle/pix
%
% Missoni Fulvio

C_tmp = F{1};    S_tmp = F{2};
[sy, sx, n_orient, n_frames, side]=size(C_tmp);
%PERMUTATION
%the convolution is computed in the {x_theta,t} domain* for each orientation
%channel (theta) and each side (left or right)
%*irrespective of y dimension. because each detectors is selective only for 
%stimulus directed ortogonaly to the its orietation 
C_tmp = permute(C_tmp,[2 4 1 3 5]);
S_tmp = permute(S_tmp,[2 4 1 3 5]);
C_tmp = reshape(C_tmp,sx,n_frames,sy*n_orient*side);
S_tmp = reshape(S_tmp,sx,n_frames,sy*n_orient*side);

if (k0*v >= 1/2) %Nyquist Limit
    error('The maximum velocity detectable by these filters is v= %d\n',num2str((1/4)/k0))
end
B=k0/3;
%sigma value
sigmat = sqrt(2*log(2)) ./ (2*pi*B);
if (n_frames <= 3*sigmat)
    error('Minimum number of frames is %d.\n',ceil(3*sigmat))
end
switch filter_type
    case 'gabor'
        filtertype = 1;
    case 'exp_decay'
        filtertype = 2;        
    case 'adelson_bergen'
        filtertype = 3;
end

n_vel=size(v,2);
%if select the adelson-bergen type then the only parameter to influence the
%velocity pref is "k"
if filtertype == 3
    kk = v;
end

%for each velocity, design the specified temporal (digital) filter
for index=1:n_vel
    if(filtertype==1) %Gabor
        %center the response of the Gabor filter in 0 to make it causal
        t2 = -floor(n_frames/2):floor(n_frames/2);
        vc=v(index);
        f0t=vc*k0;
        %relative bandwith extension
        B=k0/3;
        %sigma value
        sigmat = sqrt(2*log(2)) ./ (2*pi*B);
        pe=exp(-t2.^2/(2*sigmat^2)).*cos(2*pi*f0t.*t2);
        po=exp(-t2.^2/(2*sigmat^2)).*sin(2*pi*f0t.*t2);
        %energy normalization
        pe=pe*(1/sqrt(trapz(1,abs(pe).^2)));
        po=po*(1/sqrt(trapz(1,abs(po).^2)));
        pe(isnan(pe)) = 0;            
        pe(isinf(pe)) = 1;
        po(isnan(po)) = 0;            
        po(isinf(po)) = 1;
%         pe=pe'; po=po';
%         figure, plot(pe),hold on,plot(po)
%         pause
    end
    if filtertype==2 %Exp_Decay
        t2=0:n_frames-1;
        f0t=v(index)*k0;
        %bandwidth extension gabor-like
        B=k0/1.5;
        sigmat = sqrt(2*log(2)) ./ (2*pi*B);
        %note that sigmat is fixed in a way to obtain the same profile
        %of Gabor filters but nothing prevents to use different values
        pe=exp(-t2/(2*sigmat^2)).*cos(2*pi*f0t.*t2);
        po=exp(-t2/(2*sigmat^2)).*sin(2*pi*f0t.*t2);
        %energy normalization
        pe=pe*(1/sqrt(trapz(1,abs(pe).^2)));
        po=po*(1/sqrt(trapz(1,abs(po).^2)));
        pe(isnan(pe)) = 0;            
        pe(isinf(pe)) = 1;
        po(isnan(po)) = 0;            
        po(isinf(po)) = 1;
%         pe=pe'; po=po';
    end
    if(filtertype==3) %Adelson_Bergen
        %Biologically plausible function (BPFunctions)
        k=kk(index);
        t=0:n_frames-1; 
        n=[3 5];
        if k<0
            k=-k;
            pe = (k*t).^n(1).*exp(-k.*t).*(1/factorial(n(1))-(k.*t).^2/factorial(n(1)+2));
        else
            pe = -((k*t).^n(1)).*exp(-k.*t).*(1/factorial(n(1)) - (k.*t).^2/factorial(n(1)+2));
        end
        po=(k*t).^n(2).*exp(-k.*t).*(1/factorial(n(2))-(k.*t).^2/factorial(n(2)+2));
        %energy normalization
        pe=pe*(1/sqrt(trapz(1,abs(pe).^2)));
        po=po*(1/sqrt(trapz(1,abs(po).^2)));
        pe(isnan(pe)) = 0;            
        pe(isinf(pe)) = 1;
        po(isnan(po)) = 0;            
        po(isinf(po)) = 1;
%         pe=pe'; po=po';
    end
    %convolution
%     disp('filt_time')
% tic
    for maps=1:sy*n_orient*side
%         parfor j=1:side
            C(:,:,maps) =  (conv2(squeeze(C_tmp(:,:,maps)),pe,conv_type));
            S(:,:,maps) =  (conv2(squeeze(S_tmp(:,:,maps)),pe,conv_type));
            Ct(:,:,maps) =  (conv2(squeeze(C_tmp(:,:,maps)),po,conv_type));
            St(:,:,maps) =  (conv2(squeeze(S_tmp(:,:,maps)),po,conv_type));
%         end
    end 
% toc
    [dumb, c_frames, dumb2]=size(C);
    %resort data
    C = permute(reshape(C,sx,c_frames,sy,n_orient,side),[3 1 4 2 5]);
    S = permute(reshape(S,sx,c_frames,sy,n_orient,side),[3 1 4 2 5]);
    Ct= permute(reshape(Ct,sx,c_frames,sy,n_orient,side),[3 1 4 2 5]);
    St= permute(reshape(St,sx,c_frames,sy,n_orient,side),[3 1 4 2 5]);
    %OUTPUT (concatenate each velocity map along the sixth dimension)
    if index==1
        F_new{1} = C;
        F_new{2} = S;
        F_new{3} = Ct;
        F_new{4} = St;
    else
        F_new{1} = cat(6,F_new{1},C);
        F_new{2} = cat(6,F_new{2},S);
        F_new{3} = cat(6,F_new{3},Ct);
        F_new{4} = cat(6,F_new{4},St);
    end
    clear C S Ct St
end