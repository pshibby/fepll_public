function h = fancyfigure(varargin)
% % Function Name: fancyfigure
%
%   Same as figure but fullsize and with latex interpreter

% Citation:
% If you use this code please cite: 
% S. Parameswaran, C-A. Deledalle, L. Denis and T. Q. Nguyen, "Accelerating
% GMM-based patch priors for image restoration: Three ingredients for a 
% 100x speed-up", arXiv.
%
% License details as in license.txt
% ________________________________________

h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
set(h, 'defaulttextinterpreter', 'latex');
