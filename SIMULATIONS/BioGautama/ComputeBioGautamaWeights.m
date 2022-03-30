% close all
% clear
% clc

% load 'vel_tuning_plaid'

theta_cell_OUT = 0:pi/param.nOrient:pi-pi/param.nOrient;

[xx,tt] = meshgrid(param.prefVel,theta_cell_OUT);
tt = tt + pi*(xx<0);
xx = abs(xx);
% dx = xx.*cos(tt);
% dy = xx.*sin(tt);
[dx] = stim.vel_stim.*cos(stim.truetheta);
[dy] = stim.vel_stim.*sin(stim.truetheta);

sgmx = .5;
sgmy = .5;

CT = squeeze(e(3,:,:,:,:));
% CT = CT./max(CT,[],2);

sze = size(CT);
CT = reshape(CT,sze(1)*sze(2),[]);

dx = reshape(dx,sze(3),[]);
dy = reshape(dy,sze(3),[]);

G = [];

for k=1:size(dx,2)
    k
    for l=1:size(dx,1)
        
        g = exp(-((xx.*cos(tt)-dx(l,k)).^2/(2*sgmx^2)+...
            (xx.*sin(tt)-dy(l,k)).^2/(2*sgmy^2)));

        G = [G g(:)];
        
    end
end

CT = CT./max(CT,[],1);

F = CT*CT';
R = CT*G';
t = 100;

K = (t*eye(size(F))+F);

W = pinv(t*eye(size(F))+F)*R;
W = W';
% W = W-W.*eye(size(W));
W(W<0)=0;

% lambda = 1;
% CTnonneg=zeros(88*88);
% Gnonneg = zeros(88*88,1);
% for jj = 1:88
% for ii = 1:88
%    
%     CTnonneg((jj-1)*88+ii,[1:88]+(ii-1)*88) = CT(:,jj);
%     Gnonneg((jj-1)*88+ii,1) = G(ii,jj);
% end 
% end
% 
% WW= lsqnonneg([CTnonneg;lambda*eye(size(CTnonneg))],[Gnonneg; zeros(size(Gnonneg))]);
% WW=reshape(WW,88,[]);

supp_r = xx(:).*cos(tt(:)).*cos(tt(:)')+xx(:).*sin(tt(:)).*sin(tt(:)');
% supp_t = xx(:).*cos(tt(:)).*sin(tt(:)')-xx(:).*sin(tt(:)).*cos(tt(:)')./abs(xx(:));
supp_t = angle(exp(1i*tt(:)).*exp(-1i*tt(:)'));

W2 = exp(-(supp_r - (1.1*xx(:)')).^2/(2*0.25^2)).*exp(-(supp_t).^2/(2*(pi/2)^2));
% W2 = W2-W2.*eye(size(W2));

% tmp = load('vel_tuning_polarRDS_dur72');
% W3 = squeeze(tmp.e(3,:,:,:,:));
% W3 = reshape(W3,sze(1)*sze(2),[])';
% W3 = W3-W3.*eye(size(W3));

% W2 = (0.5+0.5*cos(2*pi/4*(xx(:).*cos(tt(:)'-tt(:)) - xx(:)')))-0*eye(size(W));
% W2 = exp(-((xx(:).*cos(tt(:)'-tt(:)) - 1.2*xx(:)').^2/(2*2^2)+(tt(:)'-tt(:)).^2/(2*(pi/2)^2))).*(0.5+0.5*cos(2*pi/4*(xx(:).*cos(tt(:)'-tt(:)) - 1.2*xx(:)')));
% W2 = exp(-((xx(:).*cos(tt(:)'-tt(:)) - 1.2*xx(:)').^2/(2*1^2)+(tt(:)'-tt(:)).^2/(2*(pi/3)^2)));

G = reshape(G,sze);

CT2 = W*CT;
CT3 = W2*CT;
CT4 = W3*CT;

CT = reshape(CT,sze);
CT2 = reshape(CT2,sze);
CT3 = reshape(CT3,sze);
CT4 = reshape(CT4,sze);
W = reshape(W,sze(1),sze(2),sze(1),sze(2));
W2 = reshape(W2,sze(1),sze(2),sze(1),sze(2));
W3 = reshape(W3,sze(1),sze(2),sze(1),sze(2));

close all
ind = 400;
% ind = 1:8;
ind_phi = 1:size(CT,4);

% % % % % % % VISUALIZATION of data
% figure,
% % % % % % % Population response at certain plaid stimulus (selected with ind)
% plot_pop_response(CT(:,:,ind,ind_phi),dx(ind,ind_phi),dy(ind,ind_phi),param.prefVel)
% figure,
% % % % % % % Desired pop response profile
% plot_pop_response(G(:,:,ind,ind_phi),dx(ind,ind_phi),dy(ind,ind_phi),param.prefVel)
% % figure,
% % % % % % % Pop Response of Higher Neural Level (Weighted with weights obtained with BioGautama on pop of plaids)
% % plot_pop_response(CT2(:,:,ind,ind_phi),dx(ind,:),dy(ind,:),param.prefVel)
% figure,
% % % % % % % Pop Response of Higher Neural Level (Weighted with imposed weights profile)
% plot_pop_response(CT3(:,:,ind,ind_phi),dx(ind,:),dy(ind,:),param.prefVel)
% figure,
% % % % % % % Pop Response of Higher Neural Level (Weighted with weights obtained with BioGautama on pop of RDS)
% plot_pop_response(CT4(:,:,ind,ind_phi),dx(ind,:),dy(ind,:),param.prefVel)
% 
% dx = xx.*cos(tt);
% dy = xx.*sin(tt);
% 
% ind = 1:8;
% figure,
% %weights obtained with BioGautama on pop of plaids
% plot_pop_response(W(:,:,ind,:),dx(ind,:),dy(ind,:),param.prefVel)
% figure,
% %weights imposed with a reasoning on intersection of constraints method
% plot_pop_response(W2(:,:,ind,:),dx(ind,:),dy(ind,:),param.prefVel)
% figure,
% %weights BioGautama on pop of RDS
% plot_pop_response(W3(:,:,ind,:),dx(ind,:),dy(ind,:),param.prefVel)
% 
% %SAVE DATA (WEIGHTS FUNCTION)
% W = reshape(W,sze(1)*sze(2),sze(1)*sze(2));
% W2 = reshape(W2,sze(1)*sze(2),sze(1)*sze(2));
% W3 = reshape(W3,sze(1)*sze(2),sze(1)*sze(2));
% 
% save("GautamaWieghts88_Plaid.mat",'W')
% save("TikhonovWieghts88.mat",'W2')
% save("GautamaWieghts88_RDS.mat",'W3')