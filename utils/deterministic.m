function varargout = deterministic(str, varargin)
%% Enter or exit deterministic code sections
%
%   state = deterministic('on')
%
%   deterministic('off', state)
%
%   state is either the current state to record or the previous
%   state to restore.


rng('default')
switch str
    case 'off'
        rng(varargin{1});
    otherwise
        varargout{1} = rng;
        rng(100);
end
