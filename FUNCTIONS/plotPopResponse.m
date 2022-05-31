
function [varargout] = plotPopResponse(response,s,varargin)

if ~isempty(varargin)&&isnumeric(varargin{1})
    rx = varargin{1};
    ry = varargin{2};
    start = 3;
else
    rx = 0;
    ry = 0;
    start = 1;
end
[theta,rho] = cart2pol(rx,ry);
m=max(rho(:));
rho = rho/(2*m);
m(isnan(rho)) = 1;
rho(isnan(rho)) = 0;

for j=1:size(rx,2)
    for i = 1:size(rx,1)
       w = .5/size(rx,2);
       h = .5/size(rx,1);
        
       posx = (1-w)*rho(i,j)*cos(theta(i,j))+0.5 - w/2;
       posy = (1-h)*rho(i,j)*sin(theta(i,j))+0.5 - h/2;
        
        %         axes('Position',[0.075+rem((i-1),8)*0.23/2 0.1+floor((i-1)/8)*0.45/2 0.2/2 0.35/2]);
        %         axes('Position',[0.075+0.5*(j-1)*0.9/9 0.1+(i-1)*0.9/8 0.9/9-0.01 0.9/8-0.01]);
        axes('Position',[posx posy w h]);
        
        
        [my_rho,my_theta] = surfPolar(response(:,:,i,j),0,90,0,s/(2*m));
        colormap(jet(256))
        %hold on maximal value
        hold on
%         [Xout, Yout] = meshgrid(x);
%         % [thetaout rhoout] = cart2pol(Xout,Yout);
%         thetaout = atan(Yout./Xout);
%         [m,indx] = max(Iout);
%         [m,indy] = max(m);
        %MAX VALUE
        [m,indOr] = max(response(:,:,i,j));
        [m,indV] = max(m);
        indOr = indOr(indV);
        vx_WTA = my_rho(indOr,indV)*cos(my_theta(indOr,indV));
        vy_WTA = my_rho(indOr,indV)*sin(my_theta(indOr,indV));
%         plot(vx_WTA,vy_WTA,'g*')
        %CENTRE OF MASS VALUE
        M = sum(sum(squeeze(response(:,:,i,j))));
        vx_CM = sum(sum(squeeze(response(:,:,i,j)).*(my_rho(1:8,:).*cos(my_theta(1:8,:))),1),2)/M;
        vy_CM = sum(sum(squeeze(response(:,:,i,j)).*(my_rho(1:8,:).*sin(my_theta(1:8,:))),1),2)/M;
%         hold on, plot(vx_CM/2,vy_CM/2,'m*')
        set(gca,'Color','none')   
        grid off
%         xlim([-3 3]);
%         ylim([-3 3]);
%         if i~=2||j~=1
%             set(gca,'xtick',[],'ytick',[]);
            axis off
%         else
%             xlabel('dx')
%             ylabel('dy')
%         end
        for k=start:numel(varargin)
            if strcmp(varargin{k},'WTA')
                varargout{k} = vx_WTA; varargout{k+1} = vy_WTA;
            end
            if strcmp(varargin{k},'CM')
                varargout{k} = vx_CM; varargout{k+1} = vy_CM;
            end
        end
    end
end