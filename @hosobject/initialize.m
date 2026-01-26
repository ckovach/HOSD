 function initialize(me,N,sampling_rate,lowpass,freqs,freqindex,varargin)

% initialize(me,buffer_size,[sampling_rate],[lowpass],[freqs],[freqindex],...,'property_name',property_value)
%
% Initialize properties for all elements in a hosobject array.
% Inputs:
%    buffer_size - window size in samples used for the HOS, filter and feature
%                  waveform estimates.
%    sampling_rate - Data sampling rate (default = 1).
%    lowpass    -  Lowpass value (default sampling_rate/2). The frequency 
%       support of HOS estimates and filters are limited below this value.
%       An appropriate lowpass setting may make a dramatic difference in
%       the efficiency of estimation, as the number of coefficients scales
%       with lowpass^me.order.
%   freqs

    if ~isscalar(N)
        X = N;
        N = size(X,1);
        me(1).do_update = true;
    else
        X = [];
     end
     if nargin >1 && ~isempty(N)
        me(1).bufferN = N;
        me(1).fftN = N;
     end
    if nargin < 6
        freqindex = [];
    end
    me(1).highpassval = 2/N;
    me(1).shighpassval = 2/N;
    if nargin > 2 && ~isempty(sampling_rate)
        me(1).sampling_rate=sampling_rate;
%                 me(1).lowpassval = me(1).lowpassval*sampling_rate;
%                 me(1).glowpassval = me(1).glowpassval*sampling_rate;
%                 me(1).highpassval = me(1).highpassval*sampling_rate;
    else
        sampling_rate = me(1).sampling_rate;
    end
    if nargin > 3 && ~isempty(lowpass)
        me(1).lowpassval=lowpass./me(1).sampling_rate;
        
        %%% By default, set the global lowpass value to the same as lowpass.
        %%% HOSD filters can only detects features whose support lies within
        %%% the filter passband, so there is normally little point in
        %%% retaining HOS coefficients involving frequencies outside the
        %%% passband.
        me(1).glowpassval=lowpass./me(1).sampling_rate;
    else
        lowpass = me(1).lowpass;
    end
     if nargin < 5 || isempty(freqs)
        freqs = me(1).fftfreq(me(1).fftN)*me(1).sampling_rate;
    end
    if isnumeric(freqs)
        freqs = {freqs};
    end
    if length(freqs) > me(1).order
        me(1).order = length(freqs);
    end

    if me(1).order > length(freqs)
        freqs(end+1:me(1).order) = freqs(end);
    end
%             me.order = order;
    me(1).freqs = freqs;
%            me(1).G = ones(sum(me(1).keepfreqs{1}),1);
    me(1).do_indexing_update = false;
    k = 1;
    me(1).do_indexing_update = false;
    while k < length(varargin)    
        me(1).(varargin{k}) = varargin{k+1};
        k=k+2;
    end
    me(1).do_indexing_update = true;
    me(1).update_frequency_indexing(freqindex,me(1).mask)
    me(1).reset();

    if length(me)>1
        me(2:end).initialize(N,sampling_rate,lowpass,freqs,freqindex,varargin{:});
    end

    if ~isempty(X)
        me.get_block(X);
    end
end