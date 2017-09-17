function gstreeplot(GStree)

p = gstree2tree(GStree, [], 0);

treeplot(p, '', 'k-')

function p = gstree2tree(GStree, p, r)

L = length(p);
p(L+1) = r;
if isfield(GStree, 'child')
    for k = 1:length(GStree.child)
        p = gstree2tree(GStree.child{k}, p, L+1);
    end
end