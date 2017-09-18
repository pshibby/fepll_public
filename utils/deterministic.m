function varargout = deterministic(str, varargin)
% % Function Name: deterministic
%
%   Enter or exit deterministic code sections
%
% Inputs:
%   str         : 'on' of 'off'
%   varargin    : state to restore when str='off'
%
% Outputs:
%   varargout   : previous state when str='on'
%
% Example:
%
%   state = deterministic('on')
%
%   %%% deterministic section here%%%
%
%   deterministic('off', state)

% Citation:
% If you use this code please cite: 
% S. Parameswaran, C-A. Deledalle, L. Denis and T. Q. Nguyen, "Accelerating
% GMM-based patch priors for image restoration: Three ingredients for a 
% 100x speed-up", arXiv.
%
% License details as in license.txt
% ________________________________________


rng('default')
switch str
    case 'off'
        rng(varargin{1});
    otherwise
        varargout{1} = rng;
        rng(100);
end
