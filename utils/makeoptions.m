function options = makeoptions(varargin)
% % Function Name: makeoptions
%
%   Extracts optional arguments
%
% Input:
%   varargin    : an even number of arguments
%                 makeoptions('arg1', value1, 'arg2', value2, ...)
%
% Output:
%   options     : a structure such that
%                 options.arg1 = value1
%                 options.arg2 = value2
%                 ...

% Citation:
% If you use this code please cite: 
% S. Parameswaran, C-A. Deledalle, L. Denis and T. Q. Nguyen, "Accelerating
% GMM-based patch priors for image restoration: Three ingredients for a 
% 100x speed-up", arXiv.
%
% License details as in license.txt
% ________________________________________



if mod(length(varargin), 2) == 1
    error('the number of options should be even');
end
options = struct();
rev = @(x) x(end:-1:1);
for k = rev(1:2:length(varargin))
    options = setfield(options, lower(varargin{floor(k/2)*2+1}), ...
                                varargin{floor(k/2)*2+2});
end
