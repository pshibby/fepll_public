function op = operators(name, M, N, varargin)
% % Function Name: operators
%
%
% Inputs:
%   name        : operator name
%                   - vignetting
%                   - id (denoising)
%                   - blur (deconvolution)
%                   - blur+border (deconv with masked borders)
%                   - motion blur
%                   - subresolution (super-resolution)
%                   - randommasking (inpainting)
%                   - cs (compressed sensing)
%
%   M           : number of rows of the image
%   N           : number of columns of the image
%   varargin    : arguments depend on the operator used (see examples)
%
% Outputs:
%   op          : stuct with function definitions and parameters:
%                   - op.isize     : size of x in y = A x
%                   - op.osize     : size of y in y = A x
%                   - op.A(x)      : apply the operator A to x
%                   - op.At(y)     : apply the adjoint A^t to y
%                   - op.AtA(x)    : apply the gram matrix A^t A to x
%                   - op.inv_AtA_plus_tauId(x, tau)
%                                  : apply (A^t A + tau Id)^-1
%                   - op.inv_AtA_plus_tauLaplacian(x, tau)
%                                  : apply (A^t A + tau Laplacian)^-1
%                   - op.Tikhonov(x, tau)
%                                  : apply (A^t A + tau Laplacian)^-1 A^t
%                   - op.A_norm2   : l2 norm || A ||_2
%                   - op.A_normF   : Frobenius norm || A ||_F
%                   - op.AtA_normF : Gram Frobenius norm || A^t A ||_F
%
% Examples:
%
%   % Gaussian blur of width 3px
%   op = operators('blur', M, N, 'width', 3')
%
%   % pre-defined blur
%   K = fspecial('motion', 10, 45);
%   op = operators('blur', M, N, 'kernel', K);
%
%   % super-resolution with 0.5px Gaussian bluring, x2 subsampling
%   % and Kaiser window appodization
%   op = operators('subresoltuion', M, N, 'width', K, 'factor', 0.5,
%                  'windowing', 'kaiser');
%
%   % inpainting with 50% of missing pixels
%   op = operators('randommasking', 'factor', 0.5);
%
%   % compressed sensing with x2 compression rate
%   op = operators('randommasking', 'factor', 0.5);

% Citation:
% If you use this code please cite:
% S. Parameswaran, C-A. Deledalle, L. Denis and T. Q. Nguyen, "Accelerating
% GMM-based patch priors for image restoration: Three ingredients for a
% 100x speed-up", arXiv.
%
% License details as in license.txt
% ________________________________________

options    = makeoptions(varargin{:});

% Laplacian

laplacian = @(x) ...
    -4 * x + ...
    x([2:end end-1], :) + ...
    x([2 1:end-1], :) + ...
    x(:, [2:end end-1]) + ...
    x(:, [2 1:end-1]);

[u v] = fftgrid(M, N);
klaplacian = zeros(M, N);
klaplacian(u == 0  & v == 0) = -4;
klaplacian(u == -1 & v == 0)  = 1;
klaplacian(u == +1 & v == 0)  = 1;
klaplacian(u == 0  & v == -1) = 1;
klaplacian(u == 0  & v == +1) = 1;
flaplacian = real(fft2(klaplacian));

op.name = name;
switch name
    case 'vignetting'
        op.isize = [M, N];
        op.osize = [M, N];

        w = 0.2+0.8*ifftshift(exp(-sqrt((u.^2/M^2 + v.^2/N^2)*6).^5));

        op.A   = @(x) w .* x;
        op.At  = @(x) w .* x;
        op.AtA = @(x) w.^2 .* x;
        op.inv_AtA_plus_tauId = @(x, tau) ...
            x ./ (w.^2 + tau);
        op.inv_AtA_plus_tauLaplacian = @(x, tau) ...
            imcgs(@(x) op.AtA(x) + tau * laplacian(x), x, ...
                  1e-6, 8, [], [], x(:) ./ w(:));
        op.Tikhonov = @(x, tau) ...
            op.inv_AtA_plus_tauLaplacian(op.At(x), tau);

        op.A_norm2   = max(w(:));
        op.A_normF   = sqrt(sum(w(:).^2));
        op.AtA_normF = sqrt(sum(w(:).^4));

    case 'id' % denoising
        op.isize = [M, N];
        op.osize = [M, N];

        op.A   = @(x) x;
        op.At  = @(x) x;
        op.AtA = @(x) x;
        op.inv_AtA_plus_tauId = @(x, tau) ...
            x / (1 + tau);
        op.inv_AtA_plus_tauLaplacian = @(x, tau) ...
            imcgs(@(x) op.AtA(x) + tau * laplacian(x), x, ...
                  1e-6, 8);
        op.Tikhonov = @(x, tau) ...
            op.inv_AtA_plus_tauLaplacian(op.At(x), tau);

        op.A_norm2   = 1;
        op.A_normF   = 1 * sqrt(M*N);
        op.AtA_normF = 1 * sqrt(M*N);

    case 'blur' % deconvolution
        K = getoptions(options, 'kernel',  []);
        if isempty(K)
            h    = getoptions(options, 'width',  1.5);
            K    = exp(-(u.^2 + v.^2) / (2 * h^2));
            MK   = ceil(8 * h) + 1;
            NK   = ceil(8 * h) + 1;
            sK    = zeros(MK, NK);
            sK(:) = ifftshift(K(-(MK-1)/2 <= u & u <= (MK-1)/2 & ...
                                -(NK-1)/2 <= v & v <= (NK-1)/2));
            K = sK / sum(sK(:));
        end
        [MK, NK] = size(K);
        if mod(MK, 2) == 0
            MK = MK + 1;
        end
        if mod(NK, 2) == 0
            NK = NK + 1;
        end
        kernel = zeros(M, N);
        kernel(fftshift(-(MK-1)/2 <= u & u <= (MK-1)/2 & ...
                        -(NK-1)/2 <= v & v <= (NK-1)/2)) = K(:);
        kernel = ifftshift(kernel);
        fkernel = fft2(kernel);

        mk = floor((MK-1)/2);
        nk = floor((NK-1)/2);

        op.taper = @(x) ...
            multiedgetaper(padarray(x, [mk nk], 'replicate', 'both'), K, 3);
        op.zp    = @(x) [zeros(mk, N) ;
                         zeros(M-2*mk, nk) x zeros(M-2*mk, nk);
                         zeros(mk, N) ];
        op.crop  = @(x) x((1+mk):(end-mk), (1+nk):(end-nk));

        op.isize = [M, N];
        op.osize = [M, N];

        op.A   = @(x) real(ifft2(fkernel .* fft2(x)));
        op.At  = @(x) real(ifft2(conj(fkernel) .* fft2(x)));
        op.AtA = @(x) real(ifft2(abs(fkernel).^2 .* fft2(x)));
        op.inv_AtA_plus_tauId = @(x, tau) ...
            real(ifft2(1./(abs(fkernel).^2 + tau) .* fft2(x)));
        op.inv_AtA_plus_tauLaplacian = @(x, tau) ...
            real(ifft2(1./(abs(fkernel).^2 - tau .* flaplacian) .* fft2(x)));
        op.Tikhonov = @(x, tau) ...
            real(ifft2(conj(fkernel)./(abs(fkernel).^2 - tau .* flaplacian) .* fft2(x)));

        op.A_norm2   = 1;
        op.A_normF   = sqrt(sum(abs(fkernel(:)).^2));
        op.AtA_normF = sqrt(sum(abs(fkernel(:)).^4));

    case 'blur+border' % deconvolution with masked borders
        K = getoptions(options, 'kernel',  [], 1);

        [MK NK] = size(K);
        kernel = zeros(M, N);
        kernel(fftshift(-(MK-1)/2 <= u & u <= (MK-1)/2 & ...
                        -(NK-1)/2 <= v & v <= (NK-1)/2)) = K;
        kernel = ifftshift(kernel);
        fkernel = fft2(kernel);

        mk = ceil((MK-1)/2);
        nk = ceil((NK-1)/2);
        op.isize = [M, N];
        op.osize = [M-2*mk, N-2*nk];

        crop = @(x) x((1+mk):(end-mk), (1+nk):(end-nk));
        zp   = @(x) [zeros(mk, N) ;
                     zeros(M-2*mk, nk) x zeros(M-2*mk, nk);
                     zeros(mk, N) ];
        op.A   = @(x) crop(ifft2(fkernel .* fft2(x)));
        op.At  = @(x) ifft2(conj(fkernel) .* fft2(zp(x)));
        op.AtA = @(x) op.At(op.A(x));
        op.inv_AtA_plus_tauId = @(x, tau) ...
            imcgs(@(x) op.AtA(x) + tau * x, x, ...
                  1e-9, 1000);
        op.inv_AtA_plus_tauLaplacian = @(x, tau) ...
            imcgs(@(x) op.AtA(x) - tau * laplacian(x), x, ...
                  1e-9, 1000, [], [], mean(x(:))*ones(size(x(:))));
        op.Tikhonov = @(x, tau) ...
            op.inv_AtA_plus_tauLaplacian(op.At(x), tau);

    case 'motionblur'
        K = fspecial('motion', 10, 45);
        op = operators('blur', M, N, 'kernel', K');
        op.name = name;
        return

    case 'subresolution' % super-resolution
        h   = getoptions(options, 'width',  0.5);
        f   = getoptions(options, 'factor',  0.7);
        windowing = getoptions(options, 'windowing', 'kaiser');

        % Define blur
        kernel = exp(-(u.^2 + v.^2) / (2 * h^2));
        kernel = kernel ./ sum(kernel(:));
        fkernel = fft2(kernel);
        op.isize = [M, N];

        % Define subsampling mask
        mask = -M/2*f <= u & u < M/2*f & ...
               -N/2*f <= v & v < N/2*f;
        op.osize = [sum(mask(:,1)), sum(mask(1,:))];
        clip = @(x) reshape(x(mask), op.osize);
        zpad = @(x) zerropadding(x, mask);

        % Define anti-aliasing windowing
        switch windowing
            case 'hamming'
                hw = mask .* max((1 - cos(2*pi*u / (M*f-1) + pi))/2 .* ...
                                 (1 - cos(2*pi*v / (N*f-1) + pi))/2, 0);
            case 'kaiser'
                alpha = -log(f);
                hw = mask.* (besseli(0, pi * alpha * sqrt(1 - (2*u/(M*f-1)).^2)) .* ...
                             besseli(0, pi * alpha * sqrt(1 - (2*v/(N*f-1)).^2)) ...
                             ./ besseli(0, pi * alpha)^2);
            otherwise
                hw = 1;
        end
        fkernel = fkernel .* hw;

        % S is the renormalization of fft2 when changing image size
        S = prod(op.osize) / prod(op.isize);

        % Main definitions
        op.A   = @(x) real(ifft2(clip(fkernel .* fft2(x))) * S);
        op.At  = @(x) real(ifft2(conj(fkernel) .* zpad(fft2(x))));
        op.AtA = @(x) real(ifft2(abs(fkernel).^2 .* mask .* fft2(x) * S));
        op.inv_AtA_plus_tauId = @(x, tau) ...
            real(ifft2(1./(mask .* abs(fkernel).^2 * S + tau) .* fft2(x)));
        op.inv_AtA_plus_tauLaplacian = @(x, tau) ...
            real(ifft2(1 ./ (mask .* abs(fkernel).^2 * S - tau * flaplacian) .* fft2(x)));
        op.Tikhonov = @(x, tau) ...
            real(ifft2(conj(fkernel) ./ (mask .* abs(fkernel).^2 * S - tau * flaplacian) ...
                       .*  zpad(fft2(x))));

        op.A_norm2   = sqrt(S);
        op.A_normF   = sqrt(sum(mask(:) .* abs(fkernel(:)).^2) * S);
        op.AtA_normF = sqrt(sum(mask(:) .* abs(fkernel(:)).^4) * S^2);

    case 'randommasking' % inpainting
        f    = getoptions(options, 'factor',  0.5);
        mask = zeros(M, N);
        idx  = randperm(M*N);
        P    = ceil(M*N*f);
        mask(idx(1:P)) = 1;
        op.isize = [M, N];
        op.osize = [M, N];
        op.mask = mask;

        op.A   = @(x) mask .* x;
        op.At  = @(x) mask .* x;
        op.AtA = @(x) mask .* x;
        op.inv_AtA_plus_tauId = @(x, tau) ...
            x ./ (mask + tau);
        op.inv_AtA_plus_tauLaplacian = @(x, tau) ...
            imcgs(@(x) op.AtA(x) + tau * laplacian(x), x, ...
                  1e-6, 8);
        op.Tikhonov = @(x, tau) ...
            op.inv_AtA_plus_tauLaplacian(op.At(x), tau);

        op.A_norm2   = 1;
        op.A_normF   = sqrt(P);
        op.AtA_normF = sqrt(P);

    case 'cs' % compressed sensing
        f    = getoptions(options, 'factor',  0.5);
        Ms   = ceil(M * sqrt(f));
        Ns   = ceil(N * sqrt(f));
        p1   = [ 1, randperm(M*N-1)+1 ];
        p2   = [ 1, randperm(M*N-1)+1 ];
        P    = Ms * Ns;
        p2   = p2(1:P);
        P    = ceil(Ms * Ns);
        f1   = @(x) reshape(x(p1), [M, N]);
        f2   = @(x) x(p2);
        f2   = @(x) f2(x(:));
        f1t  = @(x) ipa(x, M, N, p1);
        f2t  = @(x) ipa(x, M, N, p2);

        op.isize = [M, N];
        op.osize = [Ms, Ns];

        r = ones(P, 1);
        op.A   = @(x) reshape(r .* f2(idct2(f1(dct2(x)))), [Ms Ns]);
        op.At  = @(x) real(idct2(f1t(dct2(f2t(r .* reshape(x, [P 1]))))));

        ri = zeros(M, N); ri(p2) = r(:);
        op.AtA = @(x) real(idct2(f1t(dct2(ri.^2 .* idct2(f1(dct2(x)))))));

        op.inv_AtA_plus_tauId = @(x, tau) ...
            real(idct2(f1t(dct2(1./(ri.^2 + tau) .* idct2(f1(dct2(x)))))));;
        op.inv_AtA_plus_tauLaplacian = @(x, tau) ...
            imcgs(@(x) op.AtA(x) + tau * laplacian(x), x, ...
                  1e-6, 8);
        op.Tikhonov = @(x, tau) ...
            op.inv_AtA_plus_tauLaplacian(op.At(x), tau);

        op.A_norm2   = max(r);
        op.A_normF   = sqrt(mean(r.^2))*sqrt(P);
        op.AtA_normF = sqrt(mean(r.^4))*sqrt(P);

    otherwise
        error(sprintf('Operator %s not implemented', name));
end

if ~isfield(op, 'A_norm2')
    op.A_norm2 = estimate_norm2(op.A, op.At, op.isize);
end

if ~isfield(op, 'A_normF')
    op.A_normF = estimate_normF(op.A, op.isize);
end
if ~isfield(op, 'AtA_normF')
    op.AtA_normF = estimate_normF(op.AtA, op.isize);
end

function z = zerropadding(x, mask)

z = zeros(size(mask));
z(mask) = x;

function z = ipa(x, M, N, idx)

z = zeros(M*N, 1);
z(idx) = x(:);
z = reshape(z, [M, N]);

%Norm_2^2 = sqrt(max(lambdamax(A' A)))
function n2 = estimate_norm2(A, At, size)

R = 50;
vect = @(x) x(:);
z    = ones(size);
for i = 1:R
    z = At(A(z));
    z = z / norm(z);
end
n2 = sum(vect(z .* At(A(z)))) / sum(vect(z.^2));
n2 = sqrt(n2);

%Norm_F^2 = sqrt(trace(A' A))
function n2 = estimate_normF(A, size)

R = 100;
n = zeros(R, 1);
vect = @(x) x(:);
for i = 1:R
    z    = sign(randn(size));
    n(i) = norm(vect(A(z)))^2;
end
n2 = sqrt(mean(n));

function x = multiedgetaper(y, K, n)

x = y;
for i = 1:n
    x = edgetaper(x, K);
end
