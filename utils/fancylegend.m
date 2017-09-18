function h = fancylegend(varargin)
% % Function Name: fancylegend
%
%   Same as legend but with latex interpreter

% Citation:
% If you use this code please cite: 
% S. Parameswaran, C-A. Deledalle, L. Denis and T. Q. Nguyen, "Accelerating
% GMM-based patch priors for image restoration: Three ingredients for a 
% 100x speed-up", arXiv.
%
% License details as in license.txt
% ________________________________________

h = legend(varargin{:});
set(h, 'interpreter', 'latex');
