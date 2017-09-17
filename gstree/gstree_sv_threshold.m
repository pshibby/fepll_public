function [GStree, GS] = gstree_sv_threshold(GStree, GS, p)

for k = 1:GS.nmodels
    S = GS.S{k};
    e = cumsum(S) / sum(S);
    idx = find(e >= p);
    t = min(idx);
    GS.S{k}(t:end) = mean(S(t:end));
    GS.Sigma{k} = GS.U{k} * diag(GS.S{k}) * GS.U{k}';
    GS.t(k) = t;
end
GStree = gstree_sv_threshold_rec(GStree, p);

function GStree = gstree_sv_threshold_rec(GStree, p)

S = GStree.S;
e = cumsum(S) / sum(S);
idx = find(e >= p);
t = min(idx);
GStree.S(t:end) = mean(GStree.S(t:end));
GStree.Sigma = GStree.U * diag(GStree.S) * GStree.U';
GStree.t = t;

if isfield(GStree, 'child')
    for k = 1:length(GStree.child)
        GStree.child{k} = gstree_sv_threshold_rec(GStree.child{k}, p);
    end
end