function [xhat, xinit] =  fepll(y, sigma, prior_model, varargin)
% % Function Name: fepll
%
%
% Inputs:
%   y           : noisy image
%   sigma       : noise std
%   prior_model : model generated using get_prior_model.m
%   varargin    : refer to retrieve arguments for a list
%
% Outputs:
%   xhat        : restored image
%   xinit       : initialization

% Citation:
% If you use this code please cite: 
% S. Parameswaran, C-A. Deledalle, L. Denis and T. Q. Nguyen, "Accelerating
% GMM-based patch priors for image restoration: Three ingredients for a 
% 100x speed-up", arXiv.
%
% License details as in license.txt
% ________________________________________

% Retrieve arguments
options    = makeoptions(varargin{:});
op         = getoptions(options, 'operator',      []);
num_iter   = getoptions(options, 'num_iter',      5);
beta_list  = getoptions(options, 'beta_list',     [1 2.^(1 + (1:num_iter))]);
Nstep      = getoptions(options, 'Nstep',         6);
randommask = getoptions(options, 'randommask',    true);
verbose    = getoptions(options, 'verbose',       true);
clip       = getoptions(options, 'clip',          true);
autonorm   = getoptions(options, 'autonorm',      false);
P          = sqrt(prior_model.GS.dim);

% Image normalization in the range [0, 1]
m = min(y(:));
M = max(y(:));
a = 5;
if m < 0 - a*sigma || M > 1 + a*sigma
    if autonorm
        warning(['Image was not in the range [0 1], ' ...
                 'it has been automatically normalized. ' ...
                 'Note that this may alter the performances. ' ...
                 'For optimal results, please normalized the input ' ...
                 'image correctly in the range [0 1].']);
        shift = m+a*sigma;
        scale = (M-m - 2*a*sigma);
        sigma = sigma / scale;
        y     = (y - shift) / scale;
    else
        warning(['Image is not in the range [0 1]. ' ...
                 'Note that this may alter the performances. ' ...
                 'For optimal results, please normalized the input ' ...
                 'image correctly in the range [0 1], or ' ...
                 'rerun FEPLL with option "autonorm: true"']);
    end
end

% Enter into deterministic section
state = deterministic('on');

% Check if boudaries have been cropped for deblurring
has_been_croped = isfield(op, 'crop') && (norm(size(y) - op.osize) ~= 0);
if has_been_croped && isfield(op, 'name') && strcmp(op.name, 'blur')
    y = op.taper(y);
end
if ~isempty(op) && norm(size(y) - op.osize) ~= 0
    error('Dimension mismatch');
end

% Compute A^t y once for all iterations
if ~isempty(op)
    Aty  = op.At(y);
end

% Initialization and betas
tstart = tic;
if isempty(op) || strcmp(op.name, 'id')
    xinit = y;
    if sigma == 0
        return;
    end
    beta_list = beta_list ./ sigma^2;
else
    factor = op.AtA_normF^2 / prod(op.isize) / op.A_norm2^2;
    if factor > 250 * sigma^2
        factor = 250 * sigma^2;
    end
    xinit = op.Tikhonov(y, sigma^2 / factor / 5);
    beta_list = beta_list * factor / sigma^2;
    if verbose
        fprintf('FEPLL: step 0 (total time: %0.2f seconds)\n', ...
                toc(tstart));
    end
end

% Core
xhat  = xinit;
[M, N] = size(xhat);
for iter = 1:num_iter
    % Update beta
    beta = beta_list(iter);

    % Extract patches
    mask   = getmask(M, N, P, Nstep, randommask);
    ztilde = getpatches(xhat, P, mask);

    % Denoise patches
    zhat   = fepll_per_patch(ztilde, 1 / beta, prior_model, varargin{:});

    % Reproject patches
    xtilde = projpatches(zhat, M, N, mask);

    % Estimate image
    if isempty(op)
        xhat = (y + beta .* sigma^2 .* xtilde) ./ ...
               (1 + beta .* sigma^2);
    else
        xhat = op.inv_AtA_plus_tauId(Aty + beta * sigma^2 .* xtilde, ...
                                     beta * sigma^2);
    end
    if verbose
        fprintf('FEPLL: step %d (total time: %0.2f seconds)\n', ...
                iter, toc(tstart));
    end
end

% Crop result to match input size
if has_been_croped
    xhat = op.crop(xhat);
end

% Clipping
%   As used in EPLL: to use in order to reproduce their results
if clip
    xhat(xhat > 1) = 1;
    xhat(xhat < 0) = 0;
end

% Rescale to original dynamic range
if exist('scale', 'var')
    xhat  = xhat  * scale + shift;
    xinit = xinit * scale + shift;
end

% Back to stochastic mode
deterministic('off', state);

return;