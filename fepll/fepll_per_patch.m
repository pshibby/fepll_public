function zhat = fepll_per_patch(ztilde, sigma2, prior_model, varargin)
% % Function Name: fepll_per_patch
%
%
% Inputs:
%   ztilde      : patches extracted from noisy image
%   sigma2      : noise variance
%   prior_model : model generated using get_prior_model.m
%   varargin    : refer to retrieve arguments for a list
%
% Outputs:
%   zhat        : restored image patches

% Citation:
% If you use this code please cite: 
% S. Parameswaran, C-A. Deledalle, L. Denis and T. Q. Nguyen, "Accelerating
% GMM-based patch priors for image restoration: Three ingredients for a 
% 100x speed-up", arXiv.
%
% License details as in license.txt
% ________________________________________

options      = makeoptions(varargin{:});
truncation   = getoptions(options, 'truncation', isfield(prior_model.GS, 't'));
truncation   = getoptions(options, 'trunc_wiener', truncation);
searchtree   = getoptions(options, 'searchtree', true);
verbose      = getoptions(options, 'verbose', true);

P2           = size(ztilde, 1);

% Remove DC component
zdc          = mean(ztilde);
ztilde       = bsxfun(@minus, ztilde, zdc);

% Gaussian selection
if searchtree
    labels = gstree_match(ztilde, prior_model.GStree, sigma2, varargin{:});
else
    labels = gs_match(ztilde, prior_model.GS, sigma2, varargin{:});
end

% Patch estimation
U          = prior_model.GS.U;
S          = prior_model.GS.S;
lab_list   = unique(labels(:))';
zhat       = zeros(size(ztilde));
for k = lab_list
    inds = labels == k;
    if truncation
        t = prior_model.GS.t(k);

        gammaj        = S{k}(1:t) ./ (S{k}(1:t) + sigma2);
        gammaP        = S{k}(P2)  ./ (S{k}(P2) + sigma2);
        zt            = ztilde(:, inds);
        ctilde        = U{k}(:, 1:t)' * zt;
        chat          = bsxfun(@times, gammaj - gammaP,  ctilde);
        zhat(:, inds) = U{k}(:, 1:t) * chat + gammaP * zt;
    else
        gammaj        = S{k} ./ (S{k} + sigma2);
        ctilde        = U{k}' * ztilde(:, inds);
        chat          = bsxfun(@times, gammaj,  ctilde);
        zhat(:, inds) = U{k} * chat;
    end

end

% Add back DC component
zhat = bsxfun(@plus, zhat, zdc);
