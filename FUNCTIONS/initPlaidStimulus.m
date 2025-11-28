function stim = initPlaidStimulus(truetheta, theta, vpld, contr, k, samples, k_gauss, alpha, varargin)
% INITPLAIDSTIMULUS
% Generates a struct array of plaid stimulus definitions
%
% OUTPUT struct fields (per stim):
%   type
%   truetheta
%   theta_g           (2×1 motion directions)
%   c                 (2×1 contrasts for gratings)
%   vpld              (scalar plaid speed)
%   vgrat             (2×1 projected grating velocities)
%   dur               
%   mode
%   disp
%   k                 (2×1 spatial frequencies)
%   apert_rad
%   k_gauss
%   alpha             (blending factor)
%   sigma_noise
%
% After this, you can call:
%       II = generate_plaid(stim(i));

    % ---------------------------------------------------------------
    % Parse optional inputs
    % ---------------------------------------------------------------
    p = inputParser;
    addParameter(p, 'disp', 0, @(x) isscalar(x) && (x == 0 || x == 1));
    addParameter(p, 'type', "plaid", @(x) ischar(x) || isstring(x));
    addParameter(p, 'sigma_noise', 0, @(x) isscalar(x));
    parse(p, varargin{:});
    opt_disp     = p.Results.disp;
    type         = p.Results.type;
    sigma_noise  = p.Results.sigma_noise;
    
    % ---------------------------------------------------------------
    % Standardize contrast input → 2-component vector per stimulus
    % ---------------------------------------------------------------
    if size(contr, 2) < 2
        % Single column or scalar: make both gratings same contrast
        mean_contrast = 0.5;
        
        % Reshape to column vector if needed
        contr = contr(:);
        
        % If scalar, make both gratings equal contrast
        if isscalar(contr)
            c = [contr, contr];
        else
            % For each value, create equal contrast pair
            c = repmat(contr, 1, 2);
        end
        
        % Clamp to valid range [0, 1]
        c = max(0, min(1, c));
    else
        % Already 2 columns - use as is
        c = contr;
    end
    
    % Ensure contrasts are in valid range
    if any(c(:) < 0) || any(c(:) > 1)
        warning('Contrast values outside [0,1] range. Clamping...');
        c = max(0, min(1, c));
    end
    
    % ---------------------------------------------------------------
    % Standardize alpha input (will be gridded like contrast was)
    % ---------------------------------------------------------------
    alpha = alpha(:);  % Ensure column vector
    
    % ---------------------------------------------------------------
    % Standardize spatial frequency → 2-component vector
    % ---------------------------------------------------------------
    if isscalar(k)
        k0 = [k, k];
    else
        k0 = k;
    end

    % ---------------------------------------------------------------
    % Build all stimulus combinations
    % truetheta   : plaid direction
    % theta(:,1)  : g1 motion direction
    % theta(:,2)  : g2 motion direction
    % vpld        : plaid motion speed
    % c(:,1/2)    : contrast pairs
    % alpha       : blending factor (replaces old contrast variation)
    % ---------------------------------------------------------------

    % Grid with alpha instead of contrast variation
    % grid #1: grating1 properties + alpha
    [~, y1, ~, c1, alpha1] = ndgrid(truetheta, theta(:,1), vpld, c(:,1), alpha);
    y1     = pagetranspose(y1);     y1     = y1(:);
    c1     = pagetranspose(c1);     c1     = c1(:);
    alpha1 = pagetranspose(alpha1); alpha1 = alpha1(:);

    % grid #2: grating2 properties + truetheta repeated
    [truetheta_grid, y2, vel_stim, c2, alpha2] = ndgrid(truetheta, theta(:,2), vpld, c(:,2), alpha);
    truetheta_vec = pagetranspose(truetheta_grid); 
    truetheta_vec = truetheta_vec(:);
    
    y2        = pagetranspose(y2);        y2        = y2(:);
    vel_stim  = pagetranspose(vel_stim);  vel_stim  = vel_stim(:);
    c2        = pagetranspose(c2);        c2        = c2(:);
    alpha2    = pagetranspose(alpha2);    alpha2    = alpha2(:);

    % ---------------------------------------------------------------
    % Allocate output struct
    % ---------------------------------------------------------------
    num_stim = length(truetheta_vec);
    stim(num_stim) = struct();

    % ---------------------------------------------------------------
    % Fill each entry
    % ---------------------------------------------------------------
    for i = 1:num_stim

        stim(i).type       = type;
        stim(i).truetheta  = truetheta_vec(i);
        stim(i).theta_g    = [y1(i), y2(i)];        % grating motion directions
        stim(i).c          = [c1(i), c2(i)];        % grating contrasts

        % plaid velocity
        stim(i).vpld       = vel_stim(i);

        % -----------------------------------------------------------
        % Project plaid velocity onto grating directions → vgrat
        % theta_g are MOTION DIRECTIONS, so subtract pi/2 to get bar orientation
        % -----------------------------------------------------------
        v_plaid = [stim(i).vpld * cos(stim(i).truetheta);
                   stim(i).vpld * sin(stim(i).truetheta)];

        vgrat = zeros(1, 2);
        for j = 1:2
            % theta_g is motion direction, bars are perpendicular (subtract pi/2)
            n_j = [cos(stim(i).theta_g(j)); 
                   sin(stim(i).theta_g(j))];
            vgrat(j) = dot(v_plaid, n_j);
        end
        stim(i).vgrat = round(vgrat, 5);

        % -----------------------------------------------------------
        % Other plaid properties
        % -----------------------------------------------------------
        stim(i).dur       = 128;
        stim(i).mode      = 1;            % analytic plaid
        stim(i).disp      = opt_disp;

        stim(i).k         = k0;           % spatial freq
        stim(i).apert_rad = ceil(samples/2) + 2;

        stim(i).k_gauss   = k_gauss;
        stim(i).alpha     = alpha1(i);    % Use gridded alpha value
        stim(i).sigma_noise = sigma_noise;
    end
end