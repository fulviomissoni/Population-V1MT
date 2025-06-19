function varargout = surfPolar(I,az,el,normal,disp)
% I:        surface to plot
% (az,el):  view of the surface
% normal:   normalise I rows by the the midpoint of each row
% disp:     radial distance for the polar plot, for controlling the extent
%           of the surface plot in the radial direction

if nargin<2
    az=30;
    el=30;
    normal = true;
end

% I = matrix containig the value of each point in polar coordinates: row =
% orientations, cols = modulus
if normal == true
    I=I./abs(repmat(I(:,floor(size(I,2)/2)+1),[1 size(I,2)]));
end

% [rho theta]=meshgrid(linspace(-1,1,size(I,2)),linspace(0,pi,size(I,1)+1));
[rho, theta]=meshgrid(disp,linspace(0,pi,size(I,1)+1));

X = rho.*cos(theta); Y=rho.*sin(theta);
% [Xout, Yout]=meshgrid(linspace(min(disp),max(disp),101));
[Xout, Yout] = meshgrid( ...
    linspace(min(disp),max(disp),41), ...
    linspace(min(disp),max(disp),21));

% [thetaout rhoout] = cart2pol(Xout,Yout);
thetaout = atan(Yout./Xout);
thetaout = thetaout + (thetaout<0)*pi;
thetaout(isnan(thetaout)) = 0;
rhoout = Xout.*cos(thetaout) + Yout.*sin(thetaout);
Iout = interp2(rho,theta,[I; fliplr(I(1,:))],rhoout,thetaout);
% y = linspace(min(disp),max(disp),101);
% y = linspace(min(disp),max(disp),41);
% 
% % x = linspace(min(disp),max(disp),101);
% x = linspace(min(disp),max(disp),11);
% 
% % imagesc(y,x,Iout)
% [x,y] = meshgrid(x,y);
g = surf(Xout,Yout,Iout);
set(g,'FaceColor','interp','EdgeColor','none')
% set(g,'FaceColor',[1,1,1])
axis off
axis xy

% fvc=surf2patch(X,Y,-1+zeros(size(X)),[I; fliplr(I(1,:))]);
% 
% patch(fvc);
% shading interp; 
% 
% Plot Polar Grid
hold on
plot(X(:,[1 end])',Y(:,[1 end])','--','LineWidth',0.1,'Color',[0.5,0.5,0.5])
plot(X,Y,'--','LineWidth',0.1,'Color',[0.5,0.5,0.5])

% % xlabel('dx')
% % ylabel('dy')
% 
view(az,el);

flag=0;
varargout{1} = rho; varargout{2} = theta;
xlim([min(X,[],"all") max(X,[],"all")]);
ylim([min(X,[],"all") max(X,[],"all")]);
end
