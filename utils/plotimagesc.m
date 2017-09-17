function h = plotimagesc(img, varargin)
%% Display an image such that back corresponds to the minimum value
%  and white to the maximu one.
%
% Input/Output
%
%    img        a M x N array
%
% Optional arguments
%
%    see plotimage



h = plotimage(img, 'Adjust', 'auto', varargin{:});
if nargout == 0
    clear h;
end

