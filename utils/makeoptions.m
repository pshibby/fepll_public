function options = makeoptions(varargin)
%% Extracts optional arguments
%
% Input/Output
%
%    VARARGIN   an even number of arguments
%               makeoptions('arg1', value1, 'arg2', value2, ...)
%
%    OPTIONS    a structure such that
%               OPTIONS.arg1 = value1
%               OPTIONS.arg2 = value2
%               ...



if mod(length(varargin), 2) == 1
    error('the number of options should be even');
end
options = struct();
rev = @(x) x(end:-1:1);
for k = rev(1:2:length(varargin))
    options = setfield(options, lower(varargin{floor(k/2)*2+1}), ...
                                varargin{floor(k/2)*2+2});
end
