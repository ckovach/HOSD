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
%              the final residual's telescope exactly to the total.
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

xr = me.hos(1:max(getcomps)).xrec(X);

for k = 1:length(getcomps)

    cmpi = getcomps(k);
    tmp0 = hosf(cmpi);

    xbefore = X - squeeze(sum(xr(:, 1:cmpi-1, :), 2));   % what component cmpi was extracted from
    xafter  = X - squeeze(sum(xr(:, 1:cmpi, :), 2));     % ... minus its reconstruction
    if me.multivariate
        % The HOS of a multivariate component lives on the output of its
        % spatiotemporal filter (the signal hosf(cmpi) was fit on in run),
        % which is the multivariate object's xfilt, not the univariate
        % hosf(cmpi)'s (that one has a single channel's filter).
        xbefore = me.hos(cmpi).xfilt(xbefore);
        xafter  = me.hos(cmpi).xfilt(xafter);
    end

    % HOS on the residual after the component: a fresh batch estimate on
    % unshifted windows (initialize so nothing of the fit's accumulation
    % leaks into it).
    hafter = hosobject(tmp0);
    hafter.reset;
    hafter.get_input(xafter, [], false, true);

    if samewindows
        hbefore = hosobject(tmp0);
        hbefore.reset;
        hbefore.get_input(xbefore, [], false, true);
        B0 = hbefore.B;
        D0 = hbefore.D;
    else
        B0 = tmp0.B;
        D0 = tmp0.D;
    end

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
