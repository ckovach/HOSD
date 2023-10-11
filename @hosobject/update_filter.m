 function update_filter(me)

    %%% Update the current value of the filter according to the current
    %%% values in Bpart. This version reconstitutes the full HOS array; it should
    %%% therefore be possible to improve efficiency using the integration
    %%% matrices, Imat, etc.
    
        Bpart = zeros([size(me.freqindx.remap),size(me.B,3)]);
        for k = 1:length(me.Bpart)          
            for kk = 1:size(me.Bpart{k},3)
                bpart = Bpart(:,:,kk) + me.Bpart{k}(me.freqindx.remap+(kk-1)*size(me.Bpart{k},1)).*(me.freqindx.partialSymmetryRegions==k);
                bpart(me.freqindx.PDconj) = conj(bpart(me.freqindx.PDconj));
                Bpart(:,:,kk) = bpart;
            end
        end

       GG = Bpart.*me.H;

       GG(isnan(GG))=0;
       G = sum(GG(:,:,:),2);

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
