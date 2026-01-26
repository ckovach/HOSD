function power_iteration(me,X)

% Refine the feature and filter estimates obtained with HOSD through a power iteration method.
% The ability of HOSD to recover a moment-maximizing filter is only
% approximate, but it yields a starting estimate which can be quickly
% refined using a power iteration. This is implemented by updating the detection
% filter with the partial delay filter computed for the current feature estimate, 
% then recomputing the feature by aligning on the outputs of the new
% feature detection filter.

 piter_method = 'realignment';
%piter_method = 'power_iteration';

switch piter_method
    case 'realignment'
        %Pseudo power-iteration through realignment. 
        %Less precise but more stable.
        pdfilt =  me.partial_delay_filt(me.feature);

        me.filterfft =pdfilt;
        % delt = me.delay;

        % me.delay = delt;
    case 'power_iteration'
        
        %%% This is untested!
        pdfilt =  me.partial_delay_filt(real(ifft(conj(me.filterfft.*me.PSD(1:end-1)))),true,false,false);
        pdfilt = pdfilt./sqrt(sum(abs(pdfilt).^2));
        me.filterfft =pdfilt;
        %Actual power itertion. May be prone to numerical instability.
        
end

if me(1).adjust_lag 
%                            ffun = ifftshift(real(ifft(me(1).filterftlag)));                   
   ffun = ifftshift(real(ifft(me(1).filterftlag.*abs(me(1).waveftlag+eps))),1);   
   mph = sum(exp(-1i*2*pi*me(1).sampt(:)./me(1).fftN).*abs(ffun).^2)./sum(abs(ffun).^2);                   
   mph = mph./(abs(mph)+eps);
   if ~isnan(mph)
     me(1).lag = mph; % Circular shift to keep filter energy centered on the window
   end

end
[~,FXsh] = me.apply_filter(X,false);

me.wavefft = nanmean(FXsh,2);

   ffun = ifftshift(real(ifft(me(1).filterftlag.*abs(me(1).waveftlag+eps))),1);   
   mph = sum(exp(-1i*2*pi*me(1).sampt(:)./me(1).fftN).*abs(ffun).^2)./sum(abs(ffun).^2);                   
   mph = mph./(abs(mph)+eps);
   if ~isnan(mph)
     me(1).lag = mph; % Circular shift to keep filter energy centered on the window
   end

 [~,mxi] = max(real(ifft(me(1).filterftlag.*me(1).waveftlag+eps)));
if mxi~=1 && ~isnan(mxi) && me(1).do_filter_update
    me(1).filterfun = circshift(me(1).filterfun,ceil(-me(1).sampt(mxi)/2));
    me(1).feature= circshift(me(1).feature,floor(-me(1).sampt(mxi)/2));
end
        
