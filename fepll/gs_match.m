function labels = gs_match(y, GS, sig2, varargin)

options      = makeoptions(varargin{:});
truncation   = getoptions(options, 'truncation', isfield(GS, 't'));
truncation   = getoptions(options, 'trunc_match', truncation);

[d, n] = size(y);

sy2 = sum(y.^2, 1);

numMix = length(GS.S);

energy = zeros(numMix, size(y,2));
for k = 1:numMix
    iSPlusSig2   = 1 ./ (GS.S{k} + sig2);
    if ~truncation
        uy = GS.U{k}' * y;
        energy(k,:) = gmm_distance(uy, iSPlusSig2, ...
                                      -2*log(GS.wts(k)));
    else
        t = GS.t(k);

        % Energy for the t-1 first dimensions
        uy           = GS.U{k}(:,1:(t-1))' * y;
        energy(k, :)   = gmm_distance(uy, iSPlusSig2(1:(t-1)), ...
                                      -2 * log(GS.wts(k)));

        % Energy for the d-t+1 last dimensions
        uyc          = sy2 - sum(uy.^2, 1);
        energy(k, :) = energy(k, :) + ...
            uyc * iSPlusSig2(t) - (d - t + 1) * log(iSPlusSig2(t));
    end
end

[~, labels] = min(energy, [],1);