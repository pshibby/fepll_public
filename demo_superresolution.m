clear all
close all

addpathrec('.')
deterministic('on');

% Parameters
sig    = 2;

% Load and generate images (24, 20, 12, 34, 29, 49)
x      = double(imread('24.png'))/255;
[M, N] = size(x);
op     = operators('subresolution', M, N, 'width', 0.5, 'factor', 1/3);
sig    = sig/255;
y      = op.A(x) + sig * randn(op.osize);

% Load prior computed offline
prior_model = get_prior_model();

% Run FEPLL
tstart = tic;
xhat = fepll(y, sig, prior_model, 'operator', op);
t = toc(tstart);

% Display
fancyfigure;
subplot(2,2,1)
plotimage(x, 'range', [0 1]);
title('Reference HR image');
subplot(2,2,2)
xn = imresize(y, op.isize, 'nearest');
plotimage(xn, 'range', [0 1]);
title(sprintf('LR image PSNR %.1f SSIM %.3f', ...
              psnr(xn, x), ...
              ssim(xn, x)));
subplot(2,2,3)
xb = imresize(y, op.isize, 'bicubic');
plotimage(xb, 'range', [0 1]);
title(sprintf('Bicubic interp PSNR %.1f SSIM %.3f', ...
              psnr(xb, x), ...
              ssim(xb, x)));
subplot(2,2,4)
plotimage(xhat, 'range', [0 1]);
title(sprintf('Super-resolved image PSNR %.1f SSIM %.3f time %.2f', ...
              psnr(xhat, x), ...
              ssim(xhat, x), t));
linkaxes
