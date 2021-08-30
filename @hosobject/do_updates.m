 function do_updates(me,X,apply_window,use_shifted,initialize)

% do_updates(me,X,apply_window,use_shifted,initialize)
%
% Update bispectral, filter and feature estimates with new data.
% Updates are done for the respective estimates if the respective values of
% me.do_bsp_update, me.do_filter_update and me.do_wave_update are true AND
% me.do_update is true. 
%
% Values are updated as a weighted average of the sample
% estimates and the prior estimates according to me.hos_learning_rate,
% me.filter_adaptation_rate, me.hos_burnin, me.filter_burnin,
% me.window_number and the input sample size.
%
% Inputs:
%   X - input data as a column vector or me.buffersize x N array.
%   threshold - static cumulant threshold (default is me.thresh).
%   apply_window - If true and data has me.buffersize rows, applies window
%       specified in me.window. Default is set in
%       HOSOBJECT/APPLY_FILTER.
%   use_adaptive_threshold - If false, computes static cumulants from the
%       current sample, otherwise uses a pre-computed reference distribution.
% 
% See also HOSOBJECT/UPDATE_WAVEFORM HOSOBJECT/UPDATE_FILTER
% HOSOBJECT/UPDATE_BISPECTRUM HOSOBJECT/LEARNINGFUNCTION
%
% Copyright Christopher K. Kovach, University of Iowa 2018-2021

    if nargin < 5 || isempty(initialize)
        initialize = false;
    end
    if nargin < 4 || isempty(use_shifted)
        use_shifted = false;
    end

    if nargin < 3 || isempty(apply_window)
        apply_window = true;
    end
    X(end+1:me.fftN,:)=0;
    [Xfilt,FXsh] = me.apply_filter(X,apply_window,use_shifted);
    if ~me.do_update
       warning('Updating is currently disabled. Set do_update = true to enable.') 
       return
    end
%             X(end+1:me.fftN,:)=0;
%             [Xfilt,FXsh] = me.apply_filter(X,apply_window,use_shifted);
    getwin = me.update_criteria(Xfilt);

    if isempty(getwin)
        return
    end
    if me.do_bsp_update
       me.update_bispectrum(FXsh(:,getwin),initialize); 
    end
    if me.do_wave_update
      %  Xsh = real(ifft(FXsh(:,getwin)));

        me.update_waveform(FXsh(:,getwin),initialize); 
    end
    if me.do_filter_update
       me.update_filter; 

       lradj = me.learningfunction(me.filter_adaptation_rate,sum(getwin),me.filter_burnin);
       me.running_mean = me.running_mean*(1-lradj) + nanmean(Xfilt(:))*lradj;
       me.running_ssq = me.running_ssq*(1-lradj) + nanmean(Xfilt(:).^2)*lradj;

        me.running_var = me.running_ssq-me.running_mean.^2;

%                Xsrt = mean(sort(zscore(Xfilt(:,getwin)).^me.order),2);
%                Xsrt = mean(sort(((Xfilt(:,getwin)-me.running_mean)./sqrt(me.running_var)).^me.order),2);
%                  Xfilt = sort(Xfilt-me.running_mean)./sqrt(me.running_var);
%                Xfilt = sort((Xfilt-nanmean(Xfilt(:)))./nanstd(Xfilt(:)));
%                Xsrt = mean(sort(Xfilt(:,getwin)),2);
    end

    me.window_number = me.window_number+sum(getwin);
    me.outputbuffer = mean(Xfilt,2);
    me.shiftbuffer = real(ifft(mean(FXsh,2)));
    if me.do_CDF_update
        Xfsrt = (sort(Xfilt)-me.running_mean)./sqrt(me.running_var);
        for k = 1:me(1).threshold_order
            cdfb(:,k) = mean(Xfsrt.^k,2);
        end
        cdfb = cat(1,me.CDFbuffer,cdfb);
        [~,srti] = sort(cdfb(:,1));
        cdfb = cdfb(srti,:);
        cdfb(mod(me(1).window_number-1,me.CDFupsample+1)+1:me.CDFupsample+1:end,:) = []; %Keep size by decimating but avoid biased sampling of quantiles by cycling which quantiles are discarded

       if ~isempty(getwin) && any(getwin)
           lradj = me.learningfunction(me.filter_adaptation_rate,sum(getwin),me.filter_burnin);
           me.CDFbuffer = me.CDFbuffer*(1-lradj) + lradj*cdfb;

       end
    end
end