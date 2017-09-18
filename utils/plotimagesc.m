function h = plotimagesc(img, varargin)
% % Function Name: plotimagesc
%
%   Display an image such that black corresponds to the minimum value
%   and white to the maximum one.
%
% Inputs:
%   img         : a M x N array
%
% Outputs:
%   h           : handle on the created axes
%
% Optional arguments:
%   see plotimage

% Citation:
% If you use this code please cite:
% S. Parameswaran, C-A. Deledalle, L. Denis and T. Q. Nguyen, "Accelerating
% GMM-based patch priors for image restoration: Three ingredients for a
% 100x speed-up", arXiv.
%
% License details as in license.txt
% ________________________________________



h = plotimage(img, 'Adjust', 'auto', varargin{:});
if nargout == 0
    clear h;
end

