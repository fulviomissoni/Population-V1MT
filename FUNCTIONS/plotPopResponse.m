
function plot_pop_response(response,rx,ry,s)

[theta,ro] = cart2pol(rx,ry);
m=max(ro(:));
ro = ro/(2*m);
m(isnan(ro)) = 1;
ro(isnan(ro)) = 0;

for j=1:size(rx,2)
    for i = 1:size(rx,1)
       w = .5/size(rx,2);
       h = .5/size(rx,1);
        
       posx = (1-w)*ro(i,j)*cos(theta(i,j))+0.5 - w/2;
       posy = (1-h)*ro(i,j)*sin(theta(i,j))+0.5 - h/2;
        
        %         axes('Position',[0.075+rem((i-1),8)*0.23/2 0.1+floor((i-1)/8)*0.45/2 0.2/2 0.35/2]);
        %         axes('Position',[0.075+0.5*(j-1)*0.9/9 0.1+(i-1)*0.9/8 0.9/9-0.01 0.9/8-0.01]);
        axes('Position',[posx posy w h]);
        
        
        [Iout,x,y] = surfPolar(response(:,:,i,j),0,90,0,s/(2*m));
        colormap(jet(256))
        %hold on maximal value
        hold on
        [Xout, Yout] = meshgrid(x);
        % [thetaout rhoout] = cart2pol(Xout,Yout);
        thetaout = atan(Yout./Xout);
        [m,indx] = max(Iout);
        [m,indy] = max(m);
        scatter(y(indy),x(indx(indy)),150,'green','*')
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

    end
end