% Sort the data
theta1 = TestedGrat_1;
v1 = TestedVGrat_1;
theta2 = TestedGrat_2;
v2 = abs(TestedVGrat_2);
% v1 = sort(testedV_grat(1:3, 1, 1));
% theta2 = sort(testedTheta_grat(1:3, 1, 1, 2));
% v2 = sort(testedV_grat(1:3, 1, 1, 2));
% Create a polar plot
figure;
hold on;
% Plot arrows for the first set of data
for i = 1:length(theta1)
quiver(0, 0, v1(i) * cos(theta1(i)), v1(i) * sin(theta1(i)), ...
'Color', 'k', 'MaxHeadSize', 1, 'LineWidth', .5,'LineStyle','-.');
end
% Plot arrows for the second set of data with different colors
colors = lines(length(theta2)); % Generate distinct colors
for i = 1:length(theta2)
quiver(0, 0, v2(i) * cos(theta2(i)), v2(i) * sin(theta2(i)), ...
'Color', 'b', 'MaxHeadSize', 1, 'LineWidth', .5,'LineStyle','--');
end
testedTrueThetaPlaid = unique(stim(1).truetheta);
testedVplaid = unique(stim(1).vpld);
quiver(0, 0, testedVplaid*cos(testedTrueThetaPlaid), testedVplaid*sin(testedTrueThetaPlaid), ...
'Color', 'r', 'MaxHeadSize', 1, 'LineWidth', 1);
% Set equal aspect ratio
axis equal;
grid on;
title('Polar Plot with Arrows');
hold off;
xlim([-1.8,1.8]),ylim([-1.8,1.8])