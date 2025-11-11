function stim = initPlaidStimulus(truetheta, theta, vpld, contr, k, samples, k_gauss, varargin)
% INITPLAIDSTIMULUS
% Generates a struct array of plaid stimulus definitions
%
% OUTPUT struct fields (per stim):
%   type
%   truetheta
%   theta_g           (2×1 orientations)
%   contrast          (2×1 contrasts for gratings)
%   vpld              (scalar plaid speed)
%   vgrat             (2×1 projected grating velocities)
%   dur               
%   mode
%   disp
%   k                 (2×1 spatial frequencies)
%   apert_rad
%   k_gauss
%   alpha
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
    addParameter(p, 'alpha', 0.5, @(x) isscalar(x));
    addParameter(p, 'sigma_noise', 0, @(x) isscalar(x));
    parse(p, varargin{:});
    opt_disp     = p.Results.disp;
    type         = p.Results.type;
    alpha        = p.Results.alpha;
    sigma_noise  = p.Results.sigma_noise;
    
    % ---------------------------------------------------------------
    % Standardize contrast input → 2-component vector per stimulus
    % ---------------------------------------------------------------
    if size(contr, 2) < 2
        % Single column or scalar: expand to 2 columns
        % Make symmetric contrasts: one lower, one higher
        mean_contrast = 0.5;
        
        % Reshape to column vector if needed
        contr = contr(:);
        
        % Create two columns: [low_contrast, high_contrast]
        c = [mean_contrast - contr/2, mean_contrast + contr/2];
        
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
    % theta(:,1)  : g1 orientation
    % theta(:,2)  : g2 orientation
    % vpld        : plaid motion speed
    % c(:,1/2)    : contrast pairs
    % ---------------------------------------------------------------

    % grid #1: grating1 properties
    [~, y1, ~, c1] = ndgrid(truetheta, theta(:,1), vpld, c(:,1));
    y1  = pagetranspose(y1);   y1  = y1(:);
    c1  = pagetranspose(c1);   c1  = c1(:);

    % grid #2: grating2 properties + truetheta repeated
    [truetheta_grid, y2, vel_stim, c2] = ndgrid(truetheta, theta(:,2), vpld, c(:,2));
    truetheta_vec = pagetranspose(truetheta_grid); 
    truetheta_vec = truetheta_vec(:);

    y2  = pagetranspose(y2);        y2  = y2(:);
    vel_stim = pagetranspose(vel_stim); vel_stim = vel_stim(:);
    c2  = pagetranspose(c2);        c2  = c2(:);

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
        stim(i).theta_g    = [y1(i), y2(i)];        % grating orientations
        stim(i).c   = [c1(i), c2(i)];        % grating contrasts

        % plaid velocity
        stim(i).vpld       = vel_stim(i);

        % -----------------------------------------------------------
        % Project plaid velocity onto grating directions → vgrat
        % -----------------------------------------------------------
        v_plaid = [stim(i).vpld * cos(stim(i).truetheta);
           stim(i).vpld * sin(stim(i).truetheta)];

        vgrat = zeros(1, 2);
        for j = 1:2
            % Project onto orientation axis (perpendicular to bars)
            n_j = [cos(stim(i).theta_g(j)); sin(stim(i).theta_g(j))];
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
        stim(i).alpha     = alpha;
        stim(i).sigma_noise = sigma_noise;

        % -- NOTE --
        % make_plaid() is no longer needed;
        % this struct is already "plaid-ready" for generate_plaid().
    end
end
