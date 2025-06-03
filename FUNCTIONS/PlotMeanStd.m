function PlotMeanStd(x,y,varargin)

if isempty(gcf)
    figure,
end
if numel(varargin) >= 1
    color = varargin{1};
else
    color = 'b';
end
if size(y,2)==1||size(y,1)==1
    y_mean = y;
else
    y_mean = mean(y);
end
x_fill = [x, fliplr(x)];
curve1 = y_mean + std(y);
curve2 = y_mean - std(y);
inBetween = [curve1, fliplr(curve2)];
f = fill(x_fill, inBetween, color);
f.FaceAlpha = 0.15;
set(f,'EdgeColor','none')
hold on,
plot(x,y_mean,color)