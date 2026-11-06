function [c,wsum] = cumulant_correction_term(me,FXs,wgt)

% [c,wsum] = cumulant_correction_term(me,FXs,wgt)
%
% The part of the windowed moment estimate of the polyspectrum that is due
% to lower-order cumulants, estimated from the same windows. UPDATE_BISPECTRUM
% and HOS_REGRESS subtract it from the products when me.cumulant_correction
% lists the orders to remove.
%
% At each coefficient the HOS estimate is the weighted mean over windows of
%
%     P = X[b1] X[b2] ... X[b_{K-1}] conj(X[b_K]),   b_K = b1 + ... + b_{K-1},
%
% the K-th moment of the tapered windows' DFT coefficients. For a
% stationary process its expectation is the cumulant spectrum except where
% the signed frequencies of a subset of the factors sum to zero: there the
% moment also holds the products of the lower-order cumulants of the blocks
% of a partition of the factors. With a finite window "sums to zero" becomes
% "sums to within the taper's leakage width", because the tapered
% coefficients of nearby bins are correlated; on an order-4 diagonal slice,
% for instance, the rows near zero modulation carry the spectrum's leakage.
% This function returns the sum of those products,
%
%     c = sum over pi of prod_{A in pi} k(A),
%
% over the partitions pi of the K factors into two or more blocks whose
% sizes are all listed in me.cumulant_correction, where k(A) is the sample
% joint cumulant of the factors in block A. It is obtained from the sample
% moments
%
%     m(A) = sum_j w_j prod_{q in A} X_j^(s_q)[b_q]
%
% (w the window weights normalized to sum 1, s_q the conjugation of factor
% q) by the same recursion over the partitions of A,
%
%     k(A) = m(A) - sum over pi of prod_{B in pi} k(B),
%
% so cumulants of the orders that are not listed are taken to be zero:
%
%   2      - the second-order terms alone, with zero-mean coefficients
%            assumed: the products of the cross-moments of pairs of factors
%            (at order 4, m(12) m(34) + m(13) m(24) + m(14) m(23)), the
%            spectrum's leakage into the estimate. There are none at order 3.
%   [1 2]  - the window means as well.
%   1:K-1  - every lower order: the moment minus c is the sample cumulant.
%
% Listing order 1 makes the correction depend on the window means, which are
% not zero when the windows are aligned on a feature (event-related
% windowing): the correction then removes the feature's own mean.
%
% Inputs:
%   FXs - DFT of the tapered windows, fftN x nWindows: a cell array with one
%         entry (shared by every factor) or me.order entries (one per
%         factor), as UPDATE_BISPECTRUM forms the products from them.
%   wgt - weights of the windows in the HOS average (nWindows elements, or
%         one weight for every window). Windows with zero weight, or with a
%         non-finite coefficient at a bin the products use, are left out.
% Outputs:
%   c    - nIs x 1, on the scale of a weighted MEAN of the products: FFX - c
%          are the corrected products of the windows one by one.
%   wsum - the total weight of the windows used, so that sum(FFX.*wgt',2) -
%          wsum*c is the corrected weighted sum.
%
% The cross-moments of pairs of factors are read from the windows'
% cross-spectral matrix over the bins the coefficients use (one matrix
% product), the means from their weighted mean; blocks of three or more
% factors are formed per coefficient and window, as the products are.
%
% See also HOSOBJECT/UPDATE_BISPECTRUM HOSOBJECT/HOS_REGRESS
%
% Copyright Christopher K. Kovach, University of Iowa 2026

    K = me.order;
    orders = me.cumulant_correction;
    Is = me.freqindx.Is;
    c = zeros(size(Is,1),1);
    wsum = 0;
    parts = set_partitions(1:K,orders);
    parts = parts(cellfun(@numel,parts) > 1);
    if isempty(parts)
        return
    end
    if ~iscell(FXs)
        FXs = {FXs};
    end
    if any(cellfun(@(F) size(F,3),FXs) > 1)
        error('hosobject:cumulant_correction:channels', ...
            'The cumulant correction takes one channel per factor.')
    end

    %%% Only the bins the products use, and only the windows that enter
    %%% the average.
    u = unique(Is(:));
    if isscalar(wgt)
        wgt = wgt*ones(1,size(FXs{1},2));
    end
    keep = reshape(wgt,1,[]) ~= 0;
    for f = 1:numel(FXs)
        keep = keep & all(isfinite(FXs{f}(u,:)),1);
    end
    w = reshape(wgt(keep),1,[]);
    wsum = sum(w);
    w = w./wsum;
    X = cellfun(@(F) F(u,keep),FXs,'UniformOutput',false);
    [~,iu] = ismember(Is,u);   % row of each factor's bin in X

    R = cell(numel(X),numel(X),2);   % cross-spectral matrices, by input pair and conjugation
    kap = cell(1,2^K-1);             % block cumulants, by the bit mask of the block's factors
    for p = 1:numel(parts)
        t = 1;
        for b = 1:numel(parts{p})
            t = t.*block_cumulant(parts{p}{b});
        end
        c = c + t;
    end

    %%% (The nested functions share K, orders, X, w, iu, R and kap with
    %%% the function above; every other name in them is their own.)
    function kA = block_cumulant(A)
        key = sum(2.^(A-1));
        if isempty(kap{key})
            kA = block_moment(A);
            subs = set_partitions(A,orders);
            subs = subs(cellfun(@numel,subs) > 1);
            for sp = 1:numel(subs)
                tB = 1;
                for jb = 1:numel(subs{sp})
                    tB = tB.*block_cumulant(subs{sp}{jb});
                end
                kA = kA - tB;
            end
            kap{key} = kA;
        end
        kA = kap{key};
    end

    function m = block_moment(A)
        in = min(A,numel(X));   % the input of each factor
        if isscalar(A)
            m = X{in}*w.';
            m = m(iu(:,A));
            if A == K
                m = conj(m);
            end
        elseif numel(A) == 2
            cj = (A(2) == K) + 1;   % 2: the second factor is the conjugated one
            if isempty(R{in(1),in(2),cj})
                if cj == 2
                    R{in(1),in(2),cj} = (X{in(1)}.*w)*X{in(2)}';
                else
                    R{in(1),in(2),cj} = (X{in(1)}.*w)*X{in(2)}.';
                end
            end
            Rp = R{in(1),in(2),cj};
            m = Rp(iu(:,A(1)) + size(Rp,1)*(iu(:,A(2))-1));
        else
            P = 1;
            for q = 1:numel(A)
                Fq = X{in(q)}(iu(:,A(q)),:);
                if A(q) == K
                    Fq = conj(Fq);
                end
                P = P.*Fq;
            end
            m = P*w.';
        end
    end
end

function P = set_partitions(A,orders)
% Every partition of the elements of A into blocks whose sizes are listed in
% ORDERS, as a cell array of cell arrays of index vectors (the partition into
% one block included when numel(A) is listed).
    if isempty(A)
        P = {{}};
        return
    end
    rest = A(2:end);
    n = numel(rest);
    P = {};
    for s = 0:2^n-1
        pick = bitget(s,1:n) == 1;
        blk = [A(1),rest(pick)];
        if ismember(numel(blk),orders)
            sub = set_partitions(rest(~pick),orders);
            P = [P,cellfun(@(x) [{blk},x],sub,'UniformOutput',false)]; %#ok<AGROW>
        end
    end
end
