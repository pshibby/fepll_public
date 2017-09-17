function v = getoptions(options, name, v, mandatory)
%% Retrieve options parameter
%
%   v = getoptions(options, 'entry', v0, mandatory);
% is equivalent to the code:
%   if isfield(options, 'entry')
%       v = options.entry;
%   else
%       v = v0;
%   end



if nargin<3
    error('Not enough arguments.');
end
if nargin<4
    mandatory = 0;
end

name = lower(name);
if isfield(options, name)
    v = eval(['options.' name ';']);
else
    if mandatory
        error(['You have to provide the option "' name '".']);
    end
end