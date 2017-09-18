function prior_model = get_prior_model(p)
% % Function Name: get_prior_model
%
%
% Inputs:
%   p           : truncation parameter (0 1] (default: 0.95)
%   
% Outputs:
%   prior_model : A struct with truncated GMM and GMM-tree for use in fepll
%                 functions

% Citation:
% If you use this code please cite: 
% S. Parameswaran, C-A. Deledalle, L. Denis and T. Q. Nguyen, "Accelerating
% GMM-based patch priors for image restoration: Three ingredients for a 
% 100x speed-up", arXiv.
%
% License details as in license.txt
% ________________________________________

if ~exist('p', 'var')
    p = 0.95;
end
load('data/gmm.mat', 'GS');
load('data/gmm_tree.mat', 'GStree');
[GStree, GS]       = gstree_sv_threshold(GStree, GS, p);
prior_model.GS     = GS;
prior_model.GStree = GStree;
