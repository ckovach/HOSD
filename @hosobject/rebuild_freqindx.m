function freqindx = rebuild_freqindx(me)
% freqindx = rebuild_freqindx(me)
%
% The frequency indexing of an object whose freqindx is empty (it is not
% saved: Transient), rebuilt from the inputs recorded at the last
% update_frequency_indexing (freqindx_inputs), so that it is the indexing
% the stored HOS (B, D, G ...) was estimated with even if a band setting or
% the mask changed afterwards. Without a record (an object indexed by an
% older version, or given a ready indexing) the current settings are used,
% as update_frequency_indexing would. Warns when the recorded sampling rate
% or buffer size differ from the current ones (freq2index reads those from
% the object, so the rebuild then follows the current values) and when the
% rebuilt indexing does not match the length of the stored HOS.
%
% Copyright Christopher K. Kovach, University of Nebraska Medical Center 2026 

    inputs = me.freqindx_inputs;
    if isempty(inputs)
        [freqindx, ~, inputs] = me.compute_freqindx(me.mask);
        me.freqindx_inputs = inputs;
        return
    end
    if isfield(inputs, 'sampling_rate') && isfield(inputs, 'bufferN') && ...
            (inputs.sampling_rate ~= me.sampling_rate || inputs.bufferN ~= me.bufferN)
        warning('hosobject:rebuild_freqindx:settingsChanged', ...
            ['Rebuilding the frequency indexing: the sampling rate or buffer size changed after the HOS was estimated ', ...
             '(recorded %g Hz / %d samples, now %g Hz / %d); the stored HOS may no longer match.'], ...
            inputs.sampling_rate, inputs.bufferN, me.sampling_rate, me.bufferN);
    end
    freqindx = me.freqindx_from_inputs(inputs);
    nB = size(me.Bval, 1);   % coefficients (rows); a multichannel B has one column per channel along dim 3
    if nB > 1 && nB ~= size(freqindx.Is, 1) + 1
        warning('hosobject:rebuild_freqindx:mismatch', ...
            'The rebuilt frequency indexing has %d coefficients but the stored HOS has %d: a setting changed after the estimate.', ...
            size(freqindx.Is, 1) + 1, nB);
    end
end
