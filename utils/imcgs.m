function [c, flag] = imcgs(A, b, varargin)
% Solve A b = c with A hermitian positive definite
% with Conjugate Gradiant
% Unlike cgs, b, c does not need to be vectors

mat       = @(x) reshape(x, size(b));
vect      = @(x) x(:);
[c, flag] = cgs( @(x) vect(A(mat(x))), vect(b), varargin{:});
c         = mat(c);
