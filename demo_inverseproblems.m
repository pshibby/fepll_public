clear all
close all

addpath('.')
deterministic('on');

% Parameters
sigref = 2;

% Update prior model offline
prior_model = get_prior_model();

% Load and generate images
x = double(imread('59.png'))/255;
[M, N] = size(x);

% List of operators
tic
op{1} = operators('vignetting', M, N);
op{2} = operators('motionblur', M, N);
op{3} = operators('subresolution', M, N, 'width', 0.5, 'factor', 1/3);
op{4} = operators('randommasking', M, N, 'factor', 1/2);
op{5} = operators('blur', M, N, 'width', 1);
toc

% Run FEPLL
K = length(op);
for k = 1:K
    sig{k} = sigref/255;
    y{k}   = op{k}.A(x) + sig{k} * randn(op{k}.osize);

    tstart = tic;
    [xhat{k} xinit{k}] = fepll(y{k}, sig{k}, prior_model, ...
                               'autonorm', false, ...
                               'operator',     op{k});
    toc(tstart);
end

% Dislpay
fancyfigure;
j = 1;
for k = 1:K
    subplot(4,K,k);
    yr = imresize(real(y{k}), op{k}.isize, 'nearest');
    if ~strcmp(op{k}.name, 'randommasking')
        plotimage(yr, 'range', [0 1]);
    else
        img0 = imresize(real(y{k}), op{k}.isize, 'nearest');
        img1 = img0;
        img1(op{k}.mask == 0) = 1;
        img = cat(3, img1, img0, img0);
        plotimage(img, 'range', [0 1]);
    end
    title(sprintf('%s PSNR %.2f SSIM %.3f', op{k}.name, psnr(yr, x), ssim(yr, x)));

    subplot(4,K,1*K+k);
    plotimage(imcgs(op{k}.AtA, op{k}.At(y{k}), 1e-4, 20), 'range', [0 1]);
    title('Least-square');

    subplot(4,K,2*K+k);
    plotimage(xinit{k}, 'range', [0 1]);
    title('Tikhonov (Initialization)');

    subplot(4,K,3*K+k);
    plotimage(xhat{k}, 'range', [0 1]);
    title(sprintf('PSNR %.2f SSIM %.3f', psnr(xhat{k}, x), ssim(xhat{k}, x)));
end
linkaxes;
