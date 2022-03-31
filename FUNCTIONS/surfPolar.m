function varargout = surfPolar(I,az,el,normal,disp)

if nargin<2
    az=30;
    el=30;
    normal=1;
end

% I = matrix containig the value of each point in polar coordinates: row =
% orientations, cols = modulus
if normal==1
    I=I./abs(repmat(I(:,floor(size(I,2)/2)+1),[1 size(I,2)]));
end

% [rho theta]=meshgrid(linspace(-1,1,size(I,2)),linspace(0,pi,size(I,1)+1));
[rho theta]=meshgrid(disp,linspace(0,pi,size(I,1)+1));

X=rho.*cos(theta); Y=rho.*sin(theta);
[Xout Yout]=meshgrid(linspace(min(disp),max(disp),101));
% [thetaout rhoout] = cart2pol(Xout,Yout);
thetaout = atan(Yout./Xout);
thetaout = thetaout + (thetaout<0)*pi;
thetaout(isnan(thetaout)) = 0;
rhoout = Xout.*cos(thetaout) + Yout.*sin(thetaout);
Iout = interp2(rho,theta,[I; fliplr(I(1,:))],rhoout,thetaout);
y = linspace(min(disp),max(disp),101);
x = linspace(min(disp),max(disp),101);
imagesc(y,x,Iout)
axis xy

% fvc=surf2patch(X,Y,-1+zeros(size(X)),[I; fliplr(I(1,:))]);
% 
% patch(fvc);
% shading interp; 
% 
hold on
plot(X(:,[1 end])',Y(:,[1 end])','k--')
plot(X,Y,'k--')

% % xlabel('dx')
% % ylabel('dy')
% 
view(az,el);

flag=0;
varargout{1} = Iout; varargout{2} = x; varargout{3} = y;
end
