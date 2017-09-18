function [u, v] = fftgrid(n1, n2)
% % Function Name: fftgrid
%
%   like ndgrid except matrices u and v correspond to the 2d frequencies
%   of a n1 x n2 image in the encoding of fft2.
%
% Inputs:
%   n1          : first dimension
%   n2          : second dimension
%
% Outputs:
%   u           : frequency in the first dimension
%   v           : frequency in second dimension

% Citation:
% If you use this code please cite: 
% S. Parameswaran, C-A. Deledalle, L. Denis and T. Q. Nguyen, "Accelerating
% GMM-based patch priors for image restoration: Three ingredients for a 
% 100x speed-up", arXiv.
%
% License details as in license.txt
% ________________________________________

% Citation:
% If you use this code please cite: 
% S. Parameswaran, C-A. Deledalle, L. Denis and T. Q. Nguyen, "Accelerating
% GMM-based patch priors for image restoration: Three ingredients for a 
% 100x speed-up", arXiv.
%
% License details as in license.txt
% ________________________________________

[u, v] = ndgrid([0:floor(n1/2) -ceil(n1/2)+1:-1], ...
                [0:floor(n2/2) -ceil(n2/2)+1:-1]);
