function OUT = movingRDS_MS(samples_x,samples_y,dur,scale,velx,vely)
%In this function I used previously defined scripts for multiscale analysis
%on disparity detectors. Usually RDS are use stereoscopically so
%randomDotMS designs a pair of randomDot - I just take the first an moving
%horizontally. Note that that
[OUT] = randomDotMS(0, 0, samples_x, samples_y, scale);
OUT = OUT(:,:,1);
VY = vely;
VX = vely;
for i=2:dur
    % frame = circshift(OUT(:,:,i-1),[floor(vely) floor(velx)]);
    frame = circshift(OUT(:,:,1),[floor(VY) floor(VX)]);

    % x= velx-floor(velx);
    % y= vely-floor(vely);
    x= VX-floor(VX);
    y= VY-floor(VY);
    frame_x = circshift(frame,[0 ceil(x)]);
    frame_y = circshift(frame,[ceil(y) 0]);
    frame_xy = circshift(frame,[ceil(y) ceil(x)]);
    w11 = (1-x)*(1-y);
    w12 = (1-x)*y;
    w21 = x*(1-y);
    w22 = x*y;
    
    OUT(:,:,i) =    w11*frame + ...
                    w12*frame_y + ...
                    w21*frame_x + ....
                    w22*frame_xy;
    VX = velx + VX;
    VY = vely + VY;
    %normalise it

    OUT(:,:,i) = (OUT(:,:,i) - min(min(OUT(:,:,i))))./ ...
                (max(max(OUT(:,:,i))) - min(min(OUT(:,:,i))));
end
