function stim = initPlaidStimulus(truetheta,theta,vpld,contr,k,samples,k_gauss,varargin)
p = inputParser;

% Add the 'disp' parameter
addParameter(p, 'disp', 0, @(x) isscalar(x) && (x == 0 || x == 1));

% Add the 'type' parameter with a default value of 'plaid'
addParameter(p, 'type', "plaid", @(x) ischar(x) || isstring(x));

% Add the 'alpha' parameter with a default value of '0.5'
addParameter(p, 'alpha', 0.5, @(x) isscalar(x));

% Add the 'sigma_noise' parameter with a default value of '0'
addParameter(p, 'sigma_noise', 0, @(x) isscalar(x));

% Parse the input arguments
parse(p, varargin{:});

% Retrieve the results
disp = p.Results.disp;
type = p.Results.type;
alpha = p.Results.alpha;
sigma_noise = p.Results.sigma_noise;

% if isempty(varargin)
%     alpha = 0.5;       %alpha channel for transparency
% end

if size(contr,2)<2
    c = 1/4+[-contr/4,contr/4];
end
if isscalar(k)
    k0 = [k,k];
end

%for each combination of theta_g, truetheta and contrast I need a stim
%structure
%Info about first grating
[~, y, ~, c1] = ndgrid(truetheta, theta(:,1), vpld, c(:,1));
y = pagetranspose(y);   % pagetranspose to match meshgrid order
y = y(:);
c1 = pagetranspose(c1);
c1 = c1(:);

%Info about second grating
[truetheta, y2, vel_stim, c2] = ndgrid(truetheta, theta(:,2), vpld, c(:,2));
truetheta = pagetranspose(truetheta);
truetheta = truetheta(:); 
y2 = pagetranspose(y2);
y2 = y2(:);
vel_stim = pagetranspose(vel_stim);
vel_stim = vel_stim(:);
c2 = pagetranspose(c2);
c2 = c2(:);

num_stim = length(truetheta); % Assuming truetheta defines the number of stimuli

stim(num_stim) = struct(); % Preallocate the structure array

for i = 1:num_stim
    % stim(i).type = 'plaid';
    stim(i).type = type;
    stim(i).truetheta = truetheta(i); % true orientation
    stim(i).theta_g = theta; % true orientation
    % stim(i).vpld = vpld(i); % single numeric value for vpld
    % stim(i).theta_g = stim(i).truetheta + [y(i) y2(i)];
    stim(i).theta_g = [y(i) y2(i)];
    stim(i).contrast = [c1(i), c2(i)];
    stim(i).vel_stim = vel_stim(i);
    stim(i).vpld = stim(i).vel_stim; %useless????
    if size(stim(i).theta_g,2) == 2
        %vector projection
        stim(i).ori = [cos( ...
            wrapTo2Pi(stim(i).theta_g(1) - stim(i).truetheta)),... 
            cos( ...
            wrapTo2Pi(stim(i).theta_g(2) - stim(i).truetheta))];
        stim(i).vgrat = [stim(i).ori(1) .* stim(i).vel_stim, stim(i).ori(2) .* stim(i).vel_stim];
        stim(i).vgrat = round(stim(i).vgrat, 5, "decimals");
    else
        stim(i).vgrat = stim(i).vel_stim;
    end
    
    stim(i).dur = 43; % duration in frame
    stim(i).mode = 1;
    stim(i).disp = disp; % Use the parsed display option
    stim(i).k = k0; % single numeric value for k
    stim(i).apert_rad = ceil(samples/2) + 2;
    stim(i).k_gauss = k_gauss; % single numeric value for k_gauss
    stim(i).alpha = alpha;
    stim(i).sigma_noise = sigma_noise;
end