function varargout = popresponse_tiled(e,varargin)

sze = size(e); % Determine the size of the input data
if numel(sze) == 3
    sze(4) = 1;
end
if numel(sze) > 4
    error("Size of data can't be higher than 4")
end
p = inputParser;

addParameter(p, 'Titles', repmat({""}, sze(3), sze(4)), @(x) isstring(x) && all(size(x) == [sze(3), sze(4)])); % Ensure Titles is a cell array of strings with the correct size
% addParameter(p, 'Titles', repmat({""}, 1, sze(3)), @(x) isstring(x) && numel(x) == sze(3)); % Ensure Titles is a cell array of strings with the correct number of elements
parse(p, varargin{:}); % Parse the input arguments

titlestr = p.Results.Titles; % Retrieve the title from parsed input



t = tiledlayout(sze(3),sze(4),'TileSpacing','none');

for i = 1:sze(3)
for j = 1:sze(4)
nexttile,
surfPolar(squeeze(e(:,:,i,j)),0,90,'false',linspace(-2,2,11));
colormap(hot(256))
% axis equal

title(titlestr(i,j))
end
end
% cb = colorbar;  cb.Layout.Tile = 'east'; 

varargout{1} = t;

