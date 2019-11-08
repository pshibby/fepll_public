% % Denoising demo script as explained in:
%
% S. Parameswaran, C-A. Deledalle, L. Denis and T. Q. Nguyen, "Accelerating
% GMM-based patch priors for image restoration: Three ingredients for a
% 100x speed-up", arXiv.
%
% If you use any part of this software package please cite the above
% publication
%
%
% License details as in license.txt
% ________________________________________

clear all
close all

addpathrec('.')
deterministic('on');

% Parameters
sig    = 20;

% Load and generate images
x      = double(imread('24.png'))/255;
[M, N] = size(x);
sig    = sig / 255;
y      = x + sig * randn(M, N);

% Load prior computed offline
prior_model = get_prior_model();

% Run FEPLL
tstart = tic;
xhat = fepll(y, sig, prior_model);
toc(tstart);

% Display
fancyfigure;
subplot(1,3,1)
plotimage(x, 'range', [0 1]);
title('Reference image');
subplot(1,3,2)
plotimage(y, 'range', [0 1]);
title(sprintf('Noisy image PSNR %.2f SSIM %.3f', ...
    psnr(y, x), ...
    ssim(y, x)));
subplot(1,3,3)
plotimage(xhat, 'range', [0 1]);
title(sprintf('FEPLL PSNR %.2f SSIM %.3f', ...
    psnr(xhat, x), ...
    ssim(xhat, x)));

linkaxes
