 function update_filter(me)

    %%% Update the current value of the filter according to the current
    %%% values in Bpart. This version reconstitutes the full HOS array; it should
    %%% therefore be possible to improve efficiency using the integration
    %%% matrices, Imat, etc.
    
        Bpart = zeros(size(me.freqindx.remap));
        for k = 1:length(me.Bpart)                    
            Bpart = Bpart + me.Bpart{k}(me.freqindx.remap).*(me.freqindx.partialSymmetryRegions==k);
        end
        Bpart(me.freqindx.PDconj) = conj(Bpart(me.freqindx.PDconj));

        GG = Bpart.*me.H;

       GG(isnan(GG))=0;
       G = sum(GG(:,:),2);

        me.G = G(me.keepfreqs{1}(abs(me.freqs{1})<=me.lowpass(1)));

        if me.adjust_lag
           ffun = ifftshift(real(ifft(me.filterftlag.*abs(me.waveftlag+eps))));                   
           mph = sum(exp(-1i*2*pi*me.sampt(:)./me.fftN).*abs(ffun).^2)./sum(abs(ffun).^2);                   
           mph = mph./(abs(mph)+eps);
           if ~isnan(mph)
              me.lag = mph; % Circular shift to keep filter energy centered on the window
           end

        end

       %%% Also need to make sure that the output of the filter applied to
       %%% the feature waveform is centered with respect to the maximum!
       [~,mxi] = max(real(ifft(me.filterftlag.*me.waveftlag+eps)));
       if mxi~=1 && ~isnan(mxi)
%                     me.filterfun = circshift(me.filterfun,-me.sampt(mxi));
            me.waveform = circshift(me.waveform,me.sampt(mxi));
       end

    end