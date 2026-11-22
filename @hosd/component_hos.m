function [out, before, stopped] = component_hos(me, X, getcomps, samewindows, channel, progress)

% out = component_hos(obj, X, [getcomps], [samewindows])
% [out, before, stopped] = component_hos(obj, X, [getcomps], [samewindows], [channel], [progress])
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
%   channel  - multivariate fits only: the signal both terms are estimated
%              on. [] (default): the output of the component's
%              spatiotemporal filter, the signal the fit maximizes the
%              cumulant of. A channel index c: channel c of the residuals
%              (a column of X). A column of channel weights w (one per
%              column of X, e.g. a column of ICA_unmixing): the projection
%              of the residuals on w. A matrix of P such columns, one row
%              per column of X (eye(size(X,2)) for every channel): the P
%              projections in one call, which reconstructs each component
%              once for all of them. A channel with zero weight does not
%              enter a projection, so a missing (NaN) sample there is not a
%              missing sample of the projection. On one channel or
%              projection the residual after component k is again the
%              signal component k+1 was extracted from, so consecutive
%              components share that estimate as in a univariate fit.
%   progress - a function handle called after each component's estimates
%              (of every projection) as stop = progress(k, n): k components
%              of the n = numel(getcomps) asked for are done. When it returns
%              logical true the call stops after that component and returns
%              the components done so far; any other value goes on, 1
%              included (fprintf, say, returns the number of bytes it
%              printed). A function declared without outputs is called
%              without one; an anonymous function is asked for one and must
%              return a value. [] (default): none.
%
% Output:
%
%   out(k) - hosobject holding the component HOS in B with the parent's
%            normalization in D (so bicoh / bicohreduced give the component
%            polycoherence) and no bias term: the bias correction belongs to
%            each estimate, not to their difference.
%   before(k) - hosobject holding the first term by itself: the HOS of the
%            signal component k was extracted from, with its own
%            normalization and bias term (so bicoh gives the polycoherence
%            of that signal): the batch estimate on uniform windows, or
%            with samewindows false the HOS accumulated during the fit.
%   Both are 1 x numel(getcomps), and with a matrix of P > 1 columns as the
%   channel input P x numel(getcomps): out(p,k) and before(p,k) are what
%   out(k) and before(k) are with column p alone as the channel input. When
%   the progress function stopped the call they hold the components done,
%   the first of getcomps onward, and nothing of the others.
%   stopped - true when the progress function returned true (the outputs
%            then hold the components done up to that point: all of them
%            when it returned true after the last), else false.

% The components. The hos getter copies the whole array on every read, so it
% is read once.
hos = me.hos;
if nargin < 3 || isempty(getcomps)
    getcomps = 1:length(hos);
end
if islogical(getcomps)
    getcomps = find(getcomps);
end
if nargin < 4 || isempty(samewindows)
    samewindows = true;
end
if nargin < 5
    channel = [];
end
if nargin < 6
    progress = [];
end
if ~isempty(progress) && ~isa(progress, 'function_handle')
    error('hosd:component_hos:progress','progress must be a function handle, called as stop = progress(k, n).')
end
P = 1;   % the projections: the columns of the channel weights
if ~isempty(channel)
    if ~me.multivariate
        error('hosd:component_hos:channel','The channel input applies to a multivariate fit.')
    end
    if ~samewindows
        % the HOS accumulated during a multivariate fit is that of the
        % filter output, not of a channel
        error('hosd:component_hos:channel','The channel input needs samewindows = true.')
    end
    if ~ismatrix(X)
        error('hosd:component_hos:channel','The channel input takes a record: samples x channels.')
    end
    if isscalar(channel) && size(X,2) > 1
        channel = double(channel);
        if channel ~= round(channel) || channel < 1 || channel > size(X,2)
            error('hosd:component_hos:channel','channel must be a column index of X or one weight per column of X.')
        end
        W = zeros(size(X,2),1);
        W(channel) = 1;
    elseif isvector(channel) && numel(channel) == size(X,2)
        W = double(channel(:));   % one projection, its weights given as a row or a column
    else
        W = double(channel);      % one projection per column
    end
    if size(W,1) ~= size(X,2) || ~ismatrix(W) || ~all(any(W,1)) || ~all(isfinite(W(:)))
        error('hosd:component_hos:channel','channel must be a column index of X, or one weight per column of X in each of its columns.')
    end
    % The channels that enter each projection: with a zero weight a NaN would
    % still make the product NaN (0*NaN), and drop a window the channel or
    % the projection has in full.
    wuse = W ~= 0;
    P = size(W,2);
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

% The residuals come from two running sums, each step reconstructing one
% component: s, the sum of the reconstructions of components 1..m, adds them
% in sequence, as cumsum adds the slices of the reconstruction of the whole
% chain -- the same sums, bit for bit; d, the signal component m+1 is
% reconstructed from, subtracts them one at a time, as xrec (HOSOBJECT/XREC,
% MVHOSD/XREC) deflates its input for each next component -- the same
% reconstructions. So only arrays the size of X are held, never the
% reconstruction of every component at once (samples x components, and x
% channels for a multivariate fit), and only components 1..max(getcomps) are
% reconstructed.
m = 0;   % the components summed in s
d = X;
s = [];

stopped = false;   % the progress function stopped the call
outc = cell(P, length(getcomps));
beforec = cell(P, length(getcomps));
prev = [];   % the last residual's estimates: component, estimator settings, the hosobjects that hold them (one per projection)
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
    % multivariate fit estimated on the filter output: there each term is
    % the output of the component's own spatiotemporal filter (on one channel
    % or projection the residual is shared, as in a univariate fit). Nor
    % with a step longer than the window plus
    % one sample: the tail left after the last window can then outrun the
    % window, and get_input puts it through an ONLINE update first (learning
    % rates, a re-estimated filter and lag), which the estimate then depends
    % on.
    reuse = samewindows && (~me.multivariate || ~isempty(channel)) && ~isempty(prev) && prev.comp == cmpi-1 && ...
        round(hafter.poverlap*hafter.buffersize) <= hafter.buffersize+1 && isequal(prev.est, est);
    needbefore = samewindows && ~reuse;

    % Bring the sums to component cmpi, keeping the sum over the components
    % before it (sprev) for the signal cmpi was extracted from. Components
    % asked for out of order, or again, start the sums over: the same sums
    % again.
    if m > cmpi || (m == cmpi && needbefore && cmpi > 1)
        m = 0;
        d = X;
        s = [];
    end
    sprev = [];
    while m < cmpi
        m = m+1;
        sprev = [];   % not held while a component is reconstructed
        % component m's reconstruction, shaped as X (MVHOSD/XREC lists the
        % channels along the third dimension)
        r = reshape(hos(m).xrec(d), size(X));
        d = d - r;
        sprev = s;
        if m == 1
            s = r;
        else
            s = s + r;
        end
        r = []; %#ok<NASGU>
    end

    xafter = X - s;    % what is left after component cmpi
    if needbefore
        if cmpi > 1
            xbefore = X - sprev;   % what component cmpi was extracted from
        else
            % X less the sum over no component: zeros of the reconstruction's
            % class, which a single or GPU reconstruction passes on to X
            xbefore = X - zeros(size(s),'like',s);
        end
    end
    sprev = []; %#ok<NASGU>

    hafters = cell(1,P);
    for p = 1:P
        if p > 1
            hafter = hosobject(tmp0);
            hafter.reset;
        end
        if me.multivariate && ~isempty(channel)
            % One channel, or one projection, of the residuals.
            if needbefore
                xb = xbefore(:,wuse(:,p))*W(wuse(:,p),p);
            end
            xa = xafter(:,wuse(:,p))*W(wuse(:,p),p);
        elseif me.multivariate
            % The HOS of a multivariate component lives on the output of its
            % spatiotemporal filter (the signal hosf(cmpi) was fit on in run),
            % which is the multivariate object's xfilt, not the univariate
            % hosf(cmpi)'s (that one has a single channel's filter).
            if samewindows
                xb = hos(cmpi).xfilt(xbefore);
            end
            xa = hos(cmpi).xfilt(xafter);
        else
            if needbefore
                xb = xbefore;
            end
            xa = xafter;
        end

        hafter.get_input(xa, [], false, true);

        if ~samewindows
            hbefore = tmp0;
        elseif reuse
            hbefore = prev.h{p};
        else
            hbefore = hosobject(tmp0);
            hbefore.reset;
            hbefore.get_input(xb, [], false, true);
        end
        B0 = hbefore.B;
        D0 = hbefore.D;
        hafters{p} = hafter;
        if nargout > 1
            beforec{p,k} = hosobject(hbefore);
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

        outc{p,k} = tmp2;
    end
    prev = struct('comp', cmpi, 'est', {est}, 'h', {hafters});
    % released before the next component is reconstructed
    xafter = []; xbefore = []; xa = []; xb = []; %#ok<NASGU>

    if ~isempty(progress)
        % Only a logical true stops the call, not 1: fprintf, say, returns the
        % number of bytes it printed. A function declared without outputs is
        % called without one; an anonymous function is always asked for one.
        if nargout(progress) == 0
            progress(k, length(getcomps));
        else
            stopnow = progress(k, length(getcomps));
            if islogical(stopnow) && isscalar(stopnow) && stopnow
                stopped = true;
                outc = outc(:,1:k);
                beforec = beforec(:,1:k);
                break
            end
        end
    end
end

out = reshape([outc{:}], size(outc));
if nargout > 1
    before = reshape([beforec{:}], size(beforec));
end
