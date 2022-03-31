%   CREA UNA IMMAGINE DI RUMORE BIANCO E LA TRASLA LUNGO LE ASCISSE DI UNA
%   QUANTITA' DI dx PIXEL, E LUNGO LE ORDINATE DI dy PIXEL, IN MODO DA
%   SIMULARE DISPARITA' ORIZZONTALI E VERTICALI.


function [OUT true] = randomDotMS(dx, dy, sx, sy, scale)


dmax=ceil(max([abs(dx) abs(dy) 1])+3); %disparita' massima

% RUMORE BIANCO

I = randMS(sx+2*dmax, sy+2*dmax,scale);
OUT(:,:,1)=I(dmax+1:(sy+dmax),dmax+1:(sx+dmax));
    
[X,Y]=meshgrid(1:sx+2*dmax,1:sy+2*dmax);
XI=X+dx;
YI=Y+dy;
II=interp2(X,Y,I,XI,YI,'linear');

II=(II-min(min(II)))./(max(max(II))-min(min(II)));

OUT(:,:,2)=II(dmax+1:(sy+dmax),dmax+1:(sx+dmax));

true=ones(size(OUT));
true(:,:,1)=dx;
true(:,:,2)=dy;