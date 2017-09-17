function h = fancylegend(varargin)
%% Same as legend but with latex interpreter as default

h = legend(varargin{:});
set(h, 'interpreter', 'latex');
