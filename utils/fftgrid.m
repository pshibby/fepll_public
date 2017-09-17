function [u, v] = fftgrid(n1, n2)

[u, v] = ndgrid([0:floor(n1/2) -ceil(n1/2)+1:-1], ...
                [0:floor(n2/2) -ceil(n2/2)+1:-1]);
