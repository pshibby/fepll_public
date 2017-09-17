function tree_path = gstree_bfsearch(GStree, idx)

if isfield(GStree, 'idx')
    if (GStree.idx == idx)
        tree_path = [];
    else
        tree_path = false;
    end
else
    for k = 1:length(GStree.child)
        subtree_path = gstree_bfsearch(GStree.child{k}, idx);
        if islogical(subtree_path) && ~subtree_path
            tree_path = false;
        else
            tree_path = [k subtree_path];
            break;
        end
    end
end