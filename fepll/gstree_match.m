function labels = gstree_match(y, GStree, sig2, varargin)
% % Function Name: gstree_match
%
%
% Inputs:
%   y           : matrix containing image patches
%   GStree      : GMM-tree 
%   sig2        : noise variance
%
% Outputs:
%   labels      : index of Gaussian components each patch belongs to

% Citation:
% If you use this code please cite: 
% S. Parameswaran, C-A. Deledalle, L. Denis and T. Q. Nguyen, "Accelerating
% GMM-based patch priors for image restoration: Three ingredients for a 
% 100x speed-up", arXiv.
%
% License details as in license.txt
% ________________________________________

options      = makeoptions(varargin{:});
truncation   = getoptions(options, 'truncation', isfield(GStree, 't'));
truncation   = getoptions(options, 'trunc_match', truncation);

sy2    = sum(y.^2, 1);
labels = gstree_match_rec(y, sy2, GStree, sig2, truncation);

end

function labels = gstree_match_rec(y, sy2, GStree, sig2, truncation)

[d, n] = size(y);

if ~isfield(GStree, 'child')
    error('GStree must have children');
end

K = length(GStree.child);
energy = zeros(K, n);
for k = 1:K
    iSPlusSig2   = 1 ./ (GStree.child{k}.S + sig2);
    if ~truncation
        uy           = GStree.child{k}.U' * y;
        energy(k, :) = gmm_distance(uy, iSPlusSig2, ...
                                    -2 * log(GStree.child{k}.wts));
    else
        t = GStree.child{k}.t;

        % Energy for the t-1 first dimensions
        uy           = GStree.child{k}.U(:,1:(t-1))' * y;
        energy(k, :) = gmm_distance(uy, iSPlusSig2(1:(t-1)), ...
                                    -2 * log(GStree.child{k}.wts));

        % Energy for the d-t+1 last dimensions
        uyc          = sy2 - sum(uy.^2, 1);
        energy(k, :) = energy(k, :) + ...
            uyc * iSPlusSig2(t) - (d - t + 1) * log(iSPlusSig2(t));
    end
end
[~, childidx] = min(energy, [], 1);
labels = zeros(1, n);
for k = 1:K
    if sum(childidx == k) == 0
        continue;
    end
    if ~isfield(GStree.child{k}, 'child')
        labels(childidx == k) = GStree.child{k}.idx;
    else
        labels(childidx == k) = ...
            gstree_match_rec(y(:, childidx == k), ...
                             sy2(:, childidx == k), ...
                             GStree.child{k}, sig2, truncation);
    end
end

end
