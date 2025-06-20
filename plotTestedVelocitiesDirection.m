% Sort the data
theta1 = sort(testedTheta_grat(1:3, 1, 1, 1));
v1 = sort(testedV_grat(1:3, 1, 1));
theta2 = sort(testedTheta_grat(1:3, 1, 1, 2));
v2 = sort(testedV_grat(1:3, 1, 1, 2));
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
'Color', colors(i, :), 'MaxHeadSize', 1, 'LineWidth', .5,'LineStyle','--');
end
quiver(0, 0, cos(testedTrueThetaPlaid(1,1)), sin(testedTrueThetaPlaid(1,1)), ...
'Color', 'r', 'MaxHeadSize', 1, 'LineWidth', 1);
% Set equal aspect ratio
axis equal;
grid on;
title('Polar Plot with Arrows');
hold off;