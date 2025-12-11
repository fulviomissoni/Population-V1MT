function varargout = popresponse_tiled_countour(e, prefval, varargin)
    sze = size(e);
    if numel(sze) == 3
        sze(4) = 1;
    end
    if numel(sze) > 4
        error("Size of data can't be higher than 4")
    end
    
    p = inputParser;
    addParameter(p, 'Titles', repmat({""}, sze(3), sze(4)), @(x) iscell(x) && all(size(x) == [sze(3), sze(4)]));
    addParameter(p, 'nLevels', 20, @isnumeric);
    addParameter(p, 'Levels', [0,max(e,[],'all')]);

    parse(p, varargin{:});
    
    titlestr = p.Results.Titles;
    nLevels = p.Results.nLevels;
    % levels_new = p.Results.Levels;
    
    % t = tiledlayout("flow", 'TileSpacing', 'none');
    t = tiledlayout(sze(3),sze(4), 'TileSpacing', 'none');

    for i = 1:sze(3)
        for j = 1:sze(4)
            nexttile
            
            I = squeeze(e(:,:,i,j));
            
            % Coordinate polari: righe = orientazioni, colonne = raggio
            [rho, theta] = meshgrid(prefval, linspace(0, pi, size(I,1)+1));
            X = rho .* cos(theta); 
            Y = rho .* sin(theta);
            
            % Griglia cartesiana per interpolazione
            [Xout, Yout] = meshgrid( ...
                linspace(min(prefval), max(prefval), 41), ...
                linspace(min(prefval), max(prefval), 21));
            
            % Converti in polari
            thetaout = atan(Yout ./ Xout);
            thetaout = thetaout + (thetaout < 0) * pi;
            thetaout(isnan(thetaout)) = 0;
            rhoout = Xout .* cos(thetaout) + Yout .* sin(thetaout);
            
            % Interpola i dati
            Iout = interp2(rho, theta, [I; fliplr(I(1,:))], rhoout, thetaout);
            
            % Definisci livelli per contour
            % max_val = max(Iout(:));
            % if max_val > 0
            levels = linspace(min(Iout(:)), max(Iout(:)), 15);
            % else
            %     levels = 10;
            % end
            
            % Contour plot
            contourf(Xout, Yout, Iout, levels);
            axis equal tight off;
            colormap(jet);
            % view(0,-90)
            % Aggiungi griglia polare
            % hold on
            % plot(X(:,[1 end])', Y(:,[1 end])', '--', 'LineWidth', 0.1, 'Color', [0.5,0.5,0.5])
            % plot(X, Y, '--', 'LineWidth', 0.1, 'Color', [0.5,0.5,0.5])
            % hold off
            
            xlim([min(X(:)), max(X(:))]);
            ylim([min(X(:)), max(X(:))]);
            
            title(titlestr{i,j});
        end
    end
    
    if nargout > 0
        varargout{1} = t;
    end
end