function update_frequency_indexing(me,freqindx,mask)
      
% Update the frequency indexing for the hosobject according to the buffer
% size and parameters in freqindx struct. 
% See also FREQ2INDEX and FIND_PRINCIPAL_DOMAIN
%
% Copyright Christopher K. Kovach, University of Iowa 2018-2021
    if ~me(1).do_indexing_update
        return
    end
    if nargin < 3 || isempty(mask)
        mask = true;
    end
    % the band vectors, keepfreqs and (when no indexing is given) the
    % indexing itself come from compute_freqindx, which also sets me.freqs
    % when it is empty or does not match fftN
    if nargin < 2 || isempty(freqindx)
        [freqindx, keepfreqs, inputs] = me.compute_freqindx(mask);
        me.freqindx_inputs = inputs;           % recorded so a reload can rebuild exactly this indexing
    else
        [~, keepfreqs] = me.compute_freqindx(mask, false);
        me.freqindx_inputs = [];               % given a ready indexing: its inputs are unknown here (a copy inherits the source's record)
    end
    me.keepfreqs = keepfreqs;
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