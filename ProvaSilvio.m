% close all
clear
load('prova.mat')
x = linspace(-2,2,11);
alpha = pi/4;
v = 1*[cos(alpha);sin(alpha)];
theta = [linspace(0,pi,9)];
theta(end)=[];
vest = [cos(theta); sin(theta)]'*v;
% PR = exp(-(x-vest).^2);
PR = prova_pop_plaid;
TC = exp(-(x-x').^2);
LL = PR*TC;
LL = [LL(:,6:end);fliplr(LL(:,1:6))];

theta = [linspace(0,2*pi,17)];
theta(end)=[];

for i = 1:length(theta)
    MT(i,:) = mean(LL.*(cos(theta(i)-theta')));
end

MT = [fliplr(MT(9:16,1:end)) MT(1:8,2:end)];

figure, plotPopResponse(MT,x);

% thetaMT = [linspace(0,pi,9)];
% thetaMT(end)=[];
% for i = length(thetaMT)
%     MT(:,:,i) = LL.*cos(thetaMT(i)-theta');
% end