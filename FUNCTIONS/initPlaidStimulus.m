function stim = initPlaidStimulus(truetheta, theta, vpld, contr, k, samples, k_gauss, varargin)
    p = inputParser;
    addParameter(p, 'disp', 0, @(x) isscalar(x) && (x == 0 || x == 1));
    addParameter(p, 'type', "plaid", @(x) ischar(x) || isstring(x));
    addParameter(p, 'alpha', 0.5, @(x) isscalar(x));
    addParameter(p, 'sigma_noise', 0, @(x) isscalar(x));
    parse(p, varargin{:});
    opt_disp = p.Results.disp;
    type = p.Results.type;
    alpha = p.Results.alpha;
    sigma_noise = p.Results.sigma_noise;
    if size(contr, 2) < 2
        c = 1/4 + [-contr/4, contr/4];
    else
        c = contr;
    end
    if isscalar(k)
        k0 = [k, k];
    else
        k0 = k;
    end
    % Build parameter grids
    [~, y, ~, c1] = ndgrid(truetheta, theta(:,1), vpld, c(:,1));
    y = pagetranspose(y);
    y = y(:);
    c1 = pagetranspose(c1);
    c1 = c1(:);
    [truetheta_grid, y2, vel_stim, c2] = ndgrid(truetheta, theta(:,2), vpld, c(:,2));
    truetheta_vec = pagetranspose(truetheta_grid);
    truetheta_vec = truetheta_vec(:);
    y2 = pagetranspose(y2);
    y2 = y2(:);
    vel_stim = pagetranspose(vel_stim);
    vel_stim = vel_stim(:);
    c2 = pagetranspose(c2);
    c2 = c2(:);
    num_stim = length(truetheta_vec);
    stim(num_stim) = struct();
    for i = 1:num_stim
        stim(i).type = type;
        stim(i).truetheta = truetheta_vec(i);
        stim(i).theta_g = [y(i), y2(i)];
        stim(i).contrast = [c1(i), c2(i)];
        stim(i).vel_stim = vel_stim(i);
        stim(i).vpld = stim(i).vel_stim;
        % Project plaid velocity onto grating normals
        v_plaid = [stim(i).vel_stim * cos(stim(i).truetheta);
                   stim(i).vel_stim * sin(stim(i).truetheta)];
        vgrat = zeros(1, 2);
        for j = 1:2
            n_j = [cos(stim(i).theta_g(j)); sin(stim(i).theta_g(j))];
            vgrat(j) = dot(v_plaid, n_j);
        end
        stim(i).vgrat = round(vgrat, 5);
        % Other parameters
        stim(i).dur = 43;
        stim(i).mode = 1;
        stim(i).disp = opt_disp;
        stim(i).k = k0;
        stim(i).apert_rad = ceil(samples/2) + 2;
        stim(i).k_gauss = k_gauss;
        stim(i).alpha = alpha;
        stim(i).sigma_noise = sigma_noise;
    end
end