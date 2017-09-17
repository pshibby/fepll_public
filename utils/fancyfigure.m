function h = fancyfigure(varargin)
%% Same as figure but fullsize and with latex interpreter as default



h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
set(h, 'defaulttextinterpreter', 'latex');
