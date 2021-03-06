% % Deblurring demo script as explained in: 
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
sig    = 2;

% Load and generate image
x      = double(imread('50.png'))/255;
[M, N] = size(x);
op     = operators('blur', M, N, 'width', 3);
sig    = sig/255;
y      = op.A(x) + sig * randn(op.osize);

% Crop y and x for realistic experiment (no periodical effects)
y      = op.crop(y);
x      = op.crop(x);

% Load prior computed offline
prior_model = get_prior_model(0.99);

% Run FEPLL
tstart = tic;
xhat = fepll(y, sig, prior_model, 'operator', op);
t = toc(tstart);

% Display
fancyfigure;
subplot(1,3,1)
plotimage(x, 'range', [0 1]);
title('Reference HR image');
subplot(1,3,2)
plotimage(y, 'range', [0 1]);
title(sprintf('Blurry image PSNR %.1f SSIM %.3f', ...
              psnr(y, x), ...
              ssim(y, x)));
subplot(1,3,3)
plotimage(xhat, 'range', [0 1]);
title(sprintf('Deblurred image PSNR %.1f SSIM %.3f  time %.2f', ...
              psnr(xhat, x), ...
              ssim(xhat, x), t));
linkaxes
