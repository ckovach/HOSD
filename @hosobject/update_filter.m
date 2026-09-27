 function update_filter(me)

    %%% Update the current value of the filter according to the current
    %%% values in Bpart, through the integration matrices (see
    %%% INTEGRATION_MATRICES) rather than the full square-form arrays: for
    %%% each symmetry region k the part-k product is weighted by the
    %%% normalized, bias-corrected polycoherence and summed over all but the
    %%% first axis,
    %%%
    %%%     G = sum_k I{k}*(H.*Bpart{k}) + Iconj{k}*conj(H.*Bpart{k}),
    %%%
    %%% the same operation PARTIAL_DELAY_FILT applies to the sample products.
    %%% Channels are kept apart: G is nkept x 1 x nchannels, the layout
    %%% filterfft and set.nchannels use.
    %%%
    %%% Until 2026-09 this reconstituted the full arrays and, in doing so,
    %%% applied the conjugation to the running total (so every part but the
    %%% last was conjugated more than once) and, at orders above 3,
    %%% accumulated into a 2-D slice so that only the last part survived. In
    %%% the default get_block path that G is replaced by partial_delay_filt's
    %%% right away; get_input / do_updates used it as the filter.

        [Imats,Iconjmats] = me.integration_matrices;
        len = size(me.freqindx.Is,1);

        BC = me.B./(me.D+eps);
        bias = sqrt(me.BIASnum./(me.D.^2+eps));
        bias(isnan(bias)) = 0;
        BC = (abs(BC)-bias).*BC./(abs(BC)+eps);
        H = conj(BC./(me.D+eps));
        H = H(1:len,1,:);
        nchan = size(H,3);

        G = 0;
        for k = 1:length(me.Bpart)
            bp = me.Bpart{k};
            HF = reshape(H.*bp(1:len,1,:),len,nchan);
            HF(isnan(HF)) = 0;
            G = G + Imats{k}*HF + Iconjmats{k}*conj(HF);
        end
        G = reshape(G,[size(G,1) 1 nchan]);

        me.G = G(me.keepfreqs{1}(abs(me.freqs{1})<=me.lowpass(1)),:,:);

        if me.adjust_lag
           ffun = ifftshift(sum(abs(ifft(me.filterftlag.*abs(me.waveftlag+eps)).*permute(me.chanweight,[2 3 1])).^2,3),1);                   
           mph = sum(exp(-1i*2*pi*me.sampt(:)./me.fftN).*ffun)./sum(ffun);                   
           mph = mph./(abs(mph)+eps);
           if ~isnan(mph)
              me.lag = mph; % Circular shift to keep filter energy centered on the window
           end

        end

       %%% Also need to make sure that the output of the filter applied to
       %%% the feature waveform is centered with respect to the maximum!
       [~,mxi] = max(sum(real(ifft(me.filterftlag.*me.waveftlag+eps)).*permute(me.chanweight,[2 3 1]),3));
       if mxi~=1 && ~isnan(mxi)
%                     me.filterfun = circshift(me.filterfun,-me.sampt(mxi));
            me.waveform = circshift(me.waveform,-me.sampt(mxi));
       end

    end
