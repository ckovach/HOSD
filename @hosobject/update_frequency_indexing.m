function update_frequency_indexing(me,freqindx,mask)
      
% Update the frequency indexing for the hosobject according to the buffer
% size and parameters in freqindx struct. 
% See also FREQ2INDEX and FIND_PRINCIPAL_DOMAIN
%
% Copyright Christopher K. Kovach, University of Iowa 2018-2021
    if ~me(1).do_indexing_update
        return
    end
    order = me.order;
    freqs=me.freqs;
    if  isempty(freqs) || length(me(1).freqs{1})~=me(1).fftN
        freqs = me(1).fftfreq(me(1).fftN)*me(1).sampling_rate;
        me(1).freqs = repmat({freqs},1,me(1).order);
    end

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
    if nargin < 3 || isempty(mask)
        mask = true;
    end
    highpass = me.highpassval*me.sampling_rate;
    if length(highpass)< order
        highpass(end+1:order) = highpass(1);
    end

    keepfreqs={};
    for k = 1:length(me.freqs)                
        keepfreqs{k} =(abs(me.freqs{k})<=lowpass(k)&abs(me.freqs{k})>highpass(k));                 %#ok<*AGROW>
    %                 freqs{k} = freqs{k}(keepfreqs{k});
       % freqindex{k} = find(keepfreqs{k});
    end    
    me.keepfreqs = keepfreqs;
    %%% Initialize the indexing   
    if nargin < 2 || isempty(freqindx)
        freqindx = freq2index(me,freqs,order,lowpass,highpass,keepfreqs,me.pdonly,[],mask,xlowpass,xhighpass,slowpass,shighpass,me.diagonal_slice,me.include_signatures); %#ok<*PROPLC,*PROP>
    end

    me.freqindx  = freqindx;

    Z =zeros(size(me.freqindx.Is,1)+1,1);
    me.B = Z; 
    me.Bpart = {};
    me.Bpart(1:me.order) = {Z};
    me.D = Z;
     me.BIASnum=Z;
    me.G= ones(sum(me.keepfreqs{1}),1);
     me.reset;
end