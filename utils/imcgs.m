function [c, flag] = imcgs(A, b, varargin)
% % Function Name: imcgs
%
%   Same as cgs, but b and c does not need to be vectors

% Citation:
% If you use this code please cite: 
% S. Parameswaran, C-A. Deledalle, L. Denis and T. Q. Nguyen, "Accelerating
% GMM-based patch priors for image restoration: Three ingredients for a 
% 100x speed-up", arXiv.
%
% License details as in license.txt
% ________________________________________

mat       = @(x) reshape(x, size(b));
vect      = @(x) x(:);
[c, flag] = cgs( @(x) vect(A(mat(x))), vect(b), varargin{:});
c         = mat(c);
