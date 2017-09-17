function img = gstree_show_sample(GStree, r)

U = GStree.U;
S = GStree.S;
mu = GStree.mu;

X = bsxfun(@plus, mu, U * diag(sqrt(S))^2 * r);
X = X / mean(sqrt(S(:)).^2) * mean((sqrt(S)));

img = reshape(X(:, 1), [8 8]);
