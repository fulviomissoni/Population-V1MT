function OUT = movingRDS_MS(samples_x,samples_y,dur,scale,velx,vely)

[OUT] = randomDotMS(0, 0, samples_x, samples_y, scale);
OUT = OUT(:,:,1);
for i=2:dur
    frame = circshift(OUT(:,:,i-1),[floor(vely) floor(velx)]);
    x= velx-floor(velx);
    y= vely-floor(vely);
    frame_x = circshift(frame,[0 ceil(x)]);
    frame_y = circshift(frame,[ceil(y) 0]);
    frame_xy = circshift(frame,[ceil(y) ceil(x)]);
    w11 = (1-x)*(1-y);
    w12 = (1-x)*y;
    w21 = x*(1-y);
    w22 = x*y;
    
    OUT(:,:,i) = w11*frame+w12*frame_y+w21*frame_x+w22*frame_xy;
end
