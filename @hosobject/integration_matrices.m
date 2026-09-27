function [I,Iconj] = integration_matrices(me)

% [I,Iconj] = integration_matrices(me)
%
% Integration matrices for the partial polyspectral products. I{k}(i,j) is
% the number of square-form points with first-axis index i that map to
% reduced coefficient j, lie in symmetry region k (the point's first
% frequency is the k-th frequency of its principal-domain representative)
% and are not conjugates of the representative; Iconj{k} counts the
% conjugate ones. For a reduced-form vector v,
%
%     I{k}*v + Iconj{k}*conj(v)
%
% is the sum over all but the first axis of the square form of the part-k
% product, the operation PARTIAL_DELAY_FILT and UPDATE_FILTER need, without
% ever forming the full array. Built in one pass over freqindx.remap (see
% INTEGRATOR for the general form) and cached in me.Imats / me.Iconjmats;
% RESET, and hence any change of the frequency indexing, clears the cache.
%
% See also INTEGRATOR, PARTIAL_DELAY_FILT, UPDATE_FILTER.

% C. Kovach 2026

order = me.order;
if length(me.Imats) < order || length(me.Iconjmats) < order || isempty(me.Imats{1})
    fx = me.freqindx;
    remap = fx.remap;
    len = size(fx.Is,1);
    kept = find(remap <= len);                       % masked / absent points map to len+1
    s1 = mod(kept-1,size(remap,1))+1;                % first-axis subscript
    j = double(remap(kept));
    c = fx.PDconj(kept);
    reg = fx.partialSymmetryRegions(kept);
    I = cell(1,order);
    Iconj = cell(1,order);
    for k = 1:order
        sel = reg == k;
        I{k} = sparse(s1(sel & ~c),j(sel & ~c),1,size(remap,1),len);
        Iconj{k} = sparse(s1(sel & c),j(sel & c),1,size(remap,1),len);
        if me.use_gpu
            I{k} = gpuArray(I{k});
            Iconj{k} = gpuArray(Iconj{k});
        end
    end
    me.Imats = I;
    me.Iconjmats = Iconj;
end
I = me.Imats;
Iconj = me.Iconjmats;
