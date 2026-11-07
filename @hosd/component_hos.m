function out = component_hos(me, X, getcomps, samewindows)

% out = component_hos(obj, X, [getcomps], [samewindows])
%
% The higher-order spectrum attributable to each component alone: the HOS
% of the signal the component was extracted from (the residual after the
% components before it) MINUS the HOS of the residual after it, normalized
% like the former. Unlike partialbicoh, which is the HOS of the component's
% reconstruction by itself, this difference keeps the cross terms between
% the component and what is left, so the component spectra sum to the
% total (up to estimation error).
%
% Input:
%
%   obj      - hosd object, fitted
%   X        - Data the object was fitted on
%   getcomps - Components to return (indices or a logical mask; default all)
%   samewindows - true (default): estimate BOTH terms afresh on the same
%              uniform batch windows of X, so the difference is a
%              difference of like estimates and the component spectra plus
%              the final residual's telescope exactly to the total. The
%              residual after component k is the signal component k+1 was
%              extracted from, so when getcomps lists k+1 right after k and
%              the two are estimated with the same settings (a univariate
%              fit), that HOS is estimated once and serves as both terms.
%              false: use the HOS accumulated during the fit for the first
%              term (feature-aligned, exponentially weighted, outlier-
%              filtered windows; the original behavior). The second term is
%              a fresh batch estimate either way -- a batch B / D estimate
%              is invariant to shifting each window by an integer number of
%              samples, so aligning the residual's windows to its features
%              would change nothing.
%
% Output:
%
%   out(k) - hosobject holding the component HOS in B with the parent's
%            normalization in D (so bicoh / bicohreduced give the component
%            polycoherence) and no bias term: the bias correction belongs to
%            each estimate, not to their difference.

if nargin < 3 || isempty(getcomps)
    getcomps = 1:length(me.hos);
end
if islogical(getcomps)
    getcomps = find(getcomps);
end
if nargin < 4 || isempty(samewindows)
    samewindows = true;
end

% The per-component fit objects. hosf is the univariate fit of each
% component (the component itself, or, for a multivariate fit, a univariate
% object fit on the component's filter output); an object assembled from an
% hosobject array (no run) has no hosf, so fall back to the components.
if isempty(me.hosf)
    hosf = me.hos;
else
    hosf = me.hosf;
end

% Column m becomes the sum of the reconstructions of components 1..m, so
% each signal below is one subtraction rather than a sum over every
% component before it. cumsum adds the columns in sequence, as
% sum(xr(:,1:m,:),2) does: the same sums, bit for bit.
xr = cumsum(me.hos(1:max(getcomps)).xrec(X), 2);

prev = [];   % the last residual's estimate: component, estimator settings, B, D
for k = 1:length(getcomps)

    cmpi = getcomps(k);
    tmp0 = hosf(cmpi);

    % HOS on the residual after the component: a fresh batch estimate on
    % unshifted windows (initialize so nothing of the fit's accumulation
    % leaks into it).
    hafter = hosobject(tmp0);
    hafter.reset;

    % What that batch estimate (get_input -> chop_input; do_updates ->
    % apply_filter, update_bispectrum -> edf_outlier_mask,
    % cumulant_correction_term) reads from the reset copy to form B and D,
    % besides its input: the order, window, FFT length and step, the taper,
    % the coefficient indexing, the detection filter -- reset leaves it at
    % ones on keepfreqs{1} with lag 1, and it enters only through
    % apply_filter's NaN test and, at even orders, the sign of each window
    % -- with the channel and sample weights, the GPU switch, the
    % normalization, the EDF outlier rule, the cumulant correction, and the
    % sizes of the B and D it overwrites (all ones after reset). The rest is
    % set to constants by reset or read only after B and D are formed
    % (update_waveform, update_filter, the CDF update). Taken before the
    % estimate, which changes some of them.
    est = {hafter.order, hafter.buffersize, hafter.fftN, hafter.poverlap, hafter.win, ...
        hafter.freqindx.Is, hafter.keepfreqs{1}, hafter.G, hafter.chanweight, hafter.sampweight, ...
        hafter.use_gpu, hafter.normalization, hafter.integrated_magnitude_normalization, ...
        hafter.edf_outlier_threshold, hafter.edf_outlier_window_fraction, ...
        hafter.cumulant_correction, size(hafter.B), size(hafter.D)};
    % isequal ignores class and complexity (single(1) equals 1), yet a single
    % FFT length or channel weight makes B and D single, and reset leaves a B
    % or D on the GPU there: compare the classes and the complexity too.
    est = [est, cellfun(@class, est, 'UniformOutput', false), ...
        {cellfun(@isreal, est), class(hafter.B), class(hafter.D)}];

    % The signal component cmpi was extracted from is the residual after
    % component cmpi-1, so if that one's estimate was made with the same
    % settings it IS this one's first term, bit for bit. Not for a
    % multivariate fit: there each term is the output of the component's
    % own spatiotemporal filter. Nor with a step longer than the window plus
    % one sample: the tail left after the last window can then outrun the
    % window, and get_input puts it through an ONLINE update first (learning
    % rates, a re-estimated filter and lag), which the estimate then depends
    % on.
    reuse = samewindows && ~me.multivariate && ~isempty(prev) && prev.comp == cmpi-1 && ...
        round(hafter.poverlap*hafter.buffersize) <= hafter.buffersize+1 && isequal(prev.est, est);

    xafter = X - squeeze(xr(:, cmpi, :));    % what is left after component cmpi
    if samewindows && ~reuse
        if cmpi > 1
            xbefore = X - squeeze(xr(:, cmpi-1, :));   % what component cmpi was extracted from
        else
            % X less the sum over no component: zeros of the reconstruction's
            % class, which a single or GPU reconstruction passes on to X
            xbefore = X - squeeze(sum(xr(:, [], :), 2));
        end
    end
    if me.multivariate
        % The HOS of a multivariate component lives on the output of its
        % spatiotemporal filter (the signal hosf(cmpi) was fit on in run),
        % which is the multivariate object's xfilt, not the univariate
        % hosf(cmpi)'s (that one has a single channel's filter).
        if samewindows
            xbefore = me.hos(cmpi).xfilt(xbefore);
        end
        xafter  = me.hos(cmpi).xfilt(xafter);
    end

    hafter.get_input(xafter, [], false, true);

    if ~samewindows
        B0 = tmp0.B;
        D0 = tmp0.D;
    elseif reuse
        B0 = prev.B;
        D0 = prev.D;
    else
        hbefore = hosobject(tmp0);
        hbefore.reset;
        hbefore.get_input(xbefore, [], false, true);
        B0 = hbefore.B;
        D0 = hbefore.D;
    end
    prev = struct('comp', cmpi, 'est', {est}, 'B', hafter.B, 'D', hafter.D);

    tmp2 = hosobject(tmp0);
    tmp2.reset;
    % The difference of the two HOS estimates ...
    tmp2.B = B0 - hafter.B;
    % ... normalized like the signal the component came from ...
    tmp2.D = D0;
    % ... with no bias term (reset leaves BIASnum = 1, which get.bicoh
    % would subtract as a spurious 1 / |D|).
    tmp2.BIAS = zeros(size(tmp2.D));

    out(k) = tmp2; %#ok<AGROW>

end
