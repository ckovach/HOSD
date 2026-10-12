function freqindx = freqindx_from_inputs(me, inputs)
% freqindx = freqindx_from_inputs(me, inputs)
%
% The freq2index struct for a recorded set of indexing inputs (see
% compute_freqindx), taken from the process-wide cache
% (hosobject.freqindx_cache) when an identical set has been computed in this
% session, else computed and stored there. Objects with identical settings
% thereby share one copy of the index arrays.
%
% Copyright Christopher K. Kovach, University of Nebraska Medical Center 2026

    freqindx = hosobject.freqindx_cache('get', inputs);
    if isempty(freqindx)
        freqindx = freq2index(me, inputs.freqs, inputs.order, inputs.lowpass, inputs.highpass, inputs.keepfreqs, inputs.pdonly, [], ...
            inputs.mask, inputs.xlowpass, inputs.xhighpass, inputs.slowpass, inputs.shighpass, inputs.diagonal_slice, inputs.include_signatures); %#ok<*PROPLC,*PROP>
        hosobject.freqindx_cache('put', inputs, freqindx);
    end
end
