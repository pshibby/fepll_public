clear all
close all

addpathrec('.')
deterministic('on');

% Parameters
sig    = 2;

% Load and generate images (50, 28, 57, 5, 25, 32)
x      = double(imread('50.png'))/255;
[M, N] = size(x);
op     = operators('blur', M, N, 'width', 3);
sig    = sig/255;
y      = op.crop(op.A(x) + sig * randn(op.osize));
x      = op.crop(x);

% Load prior computed offline
prior_model = get_prior_model();

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
