function mask = getmask(M, N, P, r, userandom)

if ~exist('userandom', 'var')
    userandom = true;
end

% Half overlapping
if r == 1
    ho = 0;
else
    ho = floor((P - r)/2);
end

[X, Y] = ndgrid((1:r:M), (1:r:N));
[Ms, Ns] = size(X);

if userandom
    X = X + randi([-ho ho], [Ms Ns]);
    Y = Y + randi([-ho ho], [Ms Ns]);
end

X = mod(X - 1, M) + 1;
Y = mod(Y - 1, N) + 1;

if userandom
    X =  mod(randi(M) + X - 1, M) + 1;
    Y =  mod(randi(N) + Y - 1, N) + 1;
end

mask = zeros(M, N);
for i = 1:Ms
    for j = 1:Ns
        if X(i, j) <= M-P+1 && Y(i, j) <= N-P+1
            mask(X(i,j), Y(i,j)) = 1;
        end
    end
end
mask = mask == 1;

for i = 1:r:M-P+1
    mask(i, 1) = 1;
    mask(i, end-P+1) = 1;
end
mask(M-P+1, 1) = 1;
mask(M-P+1, end-P+1) = 1;
for j = 1:r:N-P+1
    mask(1, j) = 1;
    mask(end-P+1, j) = 1;
end
mask(1, N-P+1) = 1;
mask(end-P+1, N-P+1) = 1;
