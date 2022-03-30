function rn=randMS(sx, sy, scale)

% sx=128; sy=108;
sxold=sx;
if mod(sx,2^scale)~=0 sx=floor((sx)/2^scale)*(2^(scale+1)); end

syold=sy;
if mod(sy,2^scale)~=0 sy=floor((sy)/2^scale)*(2^(scale+1)); end

rn=rand(sy,sx);

for i=1:scale-1
    
    tmp=rand(sy/(2^i),sx/(2^i));
    for j=1:i tmp=expand(tmp); end
    
    rn=(tmp+0.7*rn)/1.7;
    
end

rn=rn(1:syold,1:sxold);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = expand(D)
sy = size(D,1);
sx = size(D,2);

[X Y] = meshgrid(1:(sx-1)/(2*sx-1):sx, ...
    1:(sy-1)/(2*sy-1):sy);

D = interp2(D,X,Y,'nearest');