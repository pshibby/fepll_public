function prior_model = get_prior_model(p)

if ~exist('p', 'var')
    p = 0.95;
end
load('data/gmm.mat', 'GS');
load('data/gmm_tree.mat', 'GStree');
[GStree, GS]       = gstree_sv_threshold(GStree, GS, p);
prior_model.GS     = GS;
prior_model.GStree = GStree;
