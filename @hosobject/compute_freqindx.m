function [freqindx, keepfreqs, inputs] = compute_freqindx(me, mask, wantIndex)
% [freqindx, keepfreqs, inputs] = compute_freqindx(me, mask, [wantIndex])
%
% The frequency indexing of the object (the freq2index struct) for its
% current order, buffer size, sampling rate and band settings, and the
% keepfreqs masks, and the struct of inputs that determine the indexing
% (recorded on the object as freqindx_inputs by update_frequency_indexing
% so that rebuild_freqindx can reproduce it). compute_freqindx(me, mask,
% false) returns only keepfreqs and inputs (freqindx = []).
%
% The indexing is a pure function of those settings, so identical settings
% return the SAME struct from a small process-wide cache
% (hosobject.freqindx_cache): objects fitted with the same parameters share
% one copy of the index arrays (hundreds of MB each at 4th order) instead
% of holding their own, and an object loaded from disk, where freqindx is
% not saved (Transient), gets its indexing back from the cache without
% recomputation when another object with the same settings has already
% built it.
%
% Sets me.freqs when it is empty or does not match fftN, as
% update_frequency_indexing always did; nothing else on the object changes.
%
% See also UPDATE_FREQUENCY_INDEXING, FREQ2INDEX, HOSOBJECT.FREQINDX_CACHE
%
% Copyright Christopher K. Kovach, University of Nebraska Medical Center 2026

    if nargin < 2 || isempty(mask)
        mask = true;
    end
    if nargin < 3 || isempty(wantIndex)
        wantIndex = true;
    end
    order = me.order;
    freqs = me.freqs;
    if isempty(freqs) || length(me(1).freqs{1}) ~= me(1).fftN
        freqs = me(1).fftfreq(me(1).fftN)*me(1).sampling_rate;
        me(1).freqs = repmat({freqs},1,me(1).order);
    end
    freqs = me.freqs;

    lowpass = me.lowpassval*me.sampling_rate;
    if length(lowpass)< order-1
        lowpass(end+1:order-1) = lowpass(end);
    end
    if length(lowpass)< order
        lowpass(order) = me.glowpassval*me.sampling_rate;
    end
    xlowpass = me.xlowpassval*me.sampling_rate;
    if length(xlowpass)< order
        xlowpass(end+1:order) = xlowpass(end);
    end
    slowpass = me.slowpassval*me.sampling_rate;
    shighpass = me.shighpassval*me.sampling_rate;

    xhighpass = me.xhighpassval*me.sampling_rate;
    if length(xhighpass)< order
        xhighpass(end+1:order) = xhighpass(end);
    end
    highpass = me.highpassval*me.sampling_rate;
    if length(highpass)< order
        highpass(end+1:order) = highpass(1);
    end

    keepfreqs={};
    for k = 1:length(freqs)
        keepfreqs{k} =(abs(freqs{k})<=lowpass(k)&abs(freqs{k})>highpass(k));                 %#ok<*AGROW>
    end
    % everything freq2index / find_principal_domain depend on: the explicit
    % arguments plus the sampling rate and buffer size they read from the
    % object (passtol / dgtol) and the class (subclasses may override)
    inputs = struct('class', class(me), 'order', order, 'freqs', {freqs}, 'lowpass', lowpass, 'highpass', highpass, ...
        'keepfreqs', {keepfreqs}, 'pdonly', me.pdonly, 'mask', mask, 'xlowpass', xlowpass, 'xhighpass', xhighpass, ...
        'slowpass', slowpass, 'shighpass', shighpass, 'diagonal_slice', me.diagonal_slice, 'include_signatures', me.include_signatures, ...
        'sampling_rate', me.sampling_rate, 'bufferN', me.bufferN, 'fftN', me.fftN);
    if ~wantIndex
        freqindx = [];
        return
    end
    freqindx = me.freqindx_from_inputs(inputs);
end
