function [msk,badwin,info] = edf_outlier_mask(me,A,wgt)

% [msk,badwin,info] = edf_outlier_mask(me,A,[wgt])
%
% Flag the entries of the per-window polyspectral products that dominate
% the average at their coefficient, by the effective sample size (EDF) of
% the amplitude-weighted average (Kovach 2017, "A biased look at phase
% locking", eq. 20). With the awplv normalization the polycoherence B/D at
% a coefficient is a weighted mean of unit phasors over windows t, with
% weights w_t = |X(f1) ... X(fK)| (the magnitude of the window's product),
% and
%
%     EDF = (sum_t w_t)^2 / sum_t w_t^2
%
% is its effective number of samples (the bias term sqrt(BIASnum)/D of the
% polycoherence is 1/sqrt(EDF)). Removing window t changes it by
%
%     dEDF_t = (S1 - w_t)^2 / (S2 - w_t^2) - S1^2 / S2,
%     S1 = sum_t w_t,  S2 = sum_t w_t^2,
%
% which is large and positive when window t carries a large share of the
% weight: for many windows dEDF ~ (c r)^2 - 2 c r with r = w_t / mean(w)
% and c = EDF / N, so a window whose product is ten times the mean shrinks
% the EDF by about ten effective samples whatever N is. An entry is flagged
% when dEDF exceeds me.edf_outlier_threshold; the sums are updated without
% it and the test repeats on what is left (one dominating window hides the
% next) until nothing new is flagged. A window in which more than
% me.edf_outlier_window_fraction of the coefficients is flagged is flagged
% whole (badwin, and every entry of its column in msk): an event of
% interest, sparse and strong, dominates the coefficients of its own band
% and is kept at the others; a broadband artifact dominates most of them.
%
% Inputs:
%   A   - magnitudes of the products, nIs x nWindows (x nChannels), as
%         update_bispectrum forms them (|FFX|, NaN already set to 0); an
%         entry that is 0 is never flagged.
%   wgt - optional per-window weights (nWindows x 1, or nWindows x
%         nChannels), applied before the test: the weights
%         update_bispectrum sums with.
% Outputs:
%   msk    - logical, size of A: the flagged entries.
%   badwin - logical, 1 x nWindows x nChannels: windows flagged whole.
%   info   - struct with the threshold and fraction used, the number of
%            iterations, the number of entries flagged per iteration, the
%            EDF per coefficient before and after (edf0, edf1), the
%            fraction of coefficients flagged per window (frac), the
%            windows flagged whole (badwin), the largest dEDF per window
%            before any removal (maxshrink) and whether the iteration cap
%            was reached (capped).
%
% UPDATE_BISPECTRUM applies the mask (it is the only caller), and GET_BLOCK
% drops the windows flagged whole before the alignment.
%
% See also HOSOBJECT/UPDATE_BISPECTRUM HOSOBJECT/GET_BLOCK
%
% Copyright Christopher K. Kovach, University of Iowa 2026

    maxiter = 50;
    thr = me.edf_outlier_threshold;
    wfrac = me.edf_outlier_window_fraction;
    if nargin > 2 && ~isempty(wgt) && ~isscalar(wgt)
        if isvector(wgt)
            A = A.*reshape(abs(wgt),1,[]);          % one weight per window, along dim 2
        else
            A = A.*abs(permute(wgt,[3 1 2]));       % nWindows x nChannels -> 1 x nWindows x nChannels
        end
    end
    A2 = A.^2;
    S1 = sum(A,2);
    S2 = sum(A2,2);
    info = struct('threshold',thr,'window_fraction',wfrac,'iterations',0,'nflag',zeros(1,0),'capped',false);
    info.edf0 = S1.^2./S2;
    msk = false(size(A));
    if ~(thr > 0)
        info.maxshrink = max((S1-A).^2./(S2-A2) - info.edf0,[],1);
        badwin = false([1 size(A,2) size(A,3)]);
        info.edf1 = info.edf0;
        info.frac = zeros([1 size(A,2) size(A,3)]);
        info.badwin = badwin;
        return
    end
    for it = 1:maxiter
        edf = S1.^2./S2;
        % NaN where a coefficient has one window left: not flagged. A removed
        % entry (A = 0) has a shrinkage of exactly 0, so it is never flagged
        % again.
        shrink = (S1-A).^2./(S2-A2) - edf;
        if it == 1
            info.maxshrink = max(shrink,[],1);
        end
        new = shrink > thr;
        info.nflag(it) = gather(sum(new(:)));
        info.iterations = it;
        if ~any(new(:))
            break
        end
        msk = msk | new;
        S1 = S1 - sum(A.*new,2);
        S2 = S2 - sum(A2.*new,2);
        A(new) = 0;
        A2(new) = 0;
        if it == maxiter
            info.capped = true;
            warning('hosobject:edf_outlier_mask:capped','EDF outlier flagging did not converge in %i iterations',maxiter)
        end
    end
    info.frac = mean(msk,1);
    badwin = info.frac > wfrac;
    if any(badwin(:))
        for ch = 1:size(A,3)
            msk(:,badwin(1,:,ch),ch) = true;
        end
        % the whole windows' remaining entries leave the sums too, so that
        % edf1 is the EDF of what is kept
        S1 = S1 - sum(A.*badwin,2);
        S2 = S2 - sum(A2.*badwin,2);
        info.frac = mean(msk,1);
    end
    info.edf1 = S1.^2./S2;
    info.badwin = badwin;
end
