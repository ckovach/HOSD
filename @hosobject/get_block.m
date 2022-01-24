  function  varargout = get_block(me,xin,maxiter,makeplot,segment,compno,initialize)
    
%  [A,B] = get_block(me,xin,[maxiter],[makeplot],[segment],[compno],[initialize])
%
% Applies HOSD to a block of data and applies iterative HOSD until
% convergence or maxiter is reached. This happens in the following steps:
%   1. Segment data into windows of me.buffersize length (if input is not
%      already a matrix with buffersize rows) and applies window function.
%   2. Computes polycoherence for the sample (if me.update_bispectrum = true).
%   3. Computes the partial delay filter for each window using the current
%       polycoherence estimate.
%   4. Computes the average delay filter across segments.
%   5. Aligns windows according to the delay of extrema. For even HOS orders, 
%       also entails changing sign according to the sign of the extremum.
%   6. Returns to step 4 and repeats until convergence criteria or iteration
%       limit is met.
%   7. Implements HOSD deflation by calling me(1,2:end).get_block(xresid,...) 
%      where xresid = xin - me(1,1).xrec(xin). 
% 
% Inputs:
%   xin - Input data
%   maxiter - maximum number of iterations (default = 25).
%   makeplot - If true, generate plot displaying segments and feature estimate
%              following each iteration. Default is true.
%   segment - Structure specifying how to segment the data with the
%             following fields:
%           .wint - window times
%           .Trange - time range around window times (must yield a window
%                     length of me.buffersize)
%           .fs  -  Sampling rate in which wint and Trange are specified.
%           Default values are determined from me.buffersize and me.poverlap 
%   compno - Current component number (only used for labeling in the
%            progress plot, default = 1).
%   initialize - If true, re-initialize the hosobject so that all prior
%                values are discarded. Default is true.
%
% Outputs:
%   A - A me.buffersize x N x length(me)
% See also HOSOBJECT/GET_INPUT HOSOBJECT/PARTIAL_DELAY_FILT
%          HOSOBJECT/RECONSTRUCT HOSOBJECT/XREC
%
% Copyright Christopher K. Kovach, University of Iowa 2018-2021

    if nargin < 6
        compno = 1;
    end
    Xsh=[];
    Xwin=[];

    wint=[];
    if nargin < 4 || isempty(makeplot)
        makeplot = true;
    end
    if nargin < 3 || isempty(maxiter)
        maxiter = 25;
    end
    if nargin < 5 || isempty(segment)
        segment = struct('Trange',[0 me(1).buffersize-1],'wint',[],'fs',1);
    end
    if nargin < 7 || isempty(initialize)
        initialize = true;
    end
%             if ~iscell(xin)
%                 xin = {xin};
%             end
    nxin = numel(xin);
    xisnan = isnan(xin);


    if ~isfield(segment,'fs')||isempty(segment.fs)
        segment.fs = 1;
    end

    if ~isfield(segment,'Trange')||isempty(segment.Trange)
        segment.Trange = [0 me(1).buffersize-1]./segment.fs;
        tindx = (0:me(1).buffersize-1)';
    else
         tindx = round((segment.Trange(1):1/segment.fs:segment.Trange(2))*segment.fs)';
         if length(tindx) ~= me(1).bufferN
             me(1).buffersize = length(tindx);
         end
    end


    if nxin >= me(1).bufferN
        me(1).bufferPos = 0; % Discard the buffer
        if  size(xin,1) ~=me(1).fftN %me(1).bufferN 
            if ~isfield(segment,'wint') || isempty(segment.wint)
                stepn = round(me(1).poverlap*me(1).bufferN);
                nget = nxin - me(1).bufferN+1;
                wint = (0:stepn:nget-1);%./segment.fs;
                segment.wint=wint/segment.fs;
            else
                wint = round(segment.wint*segment.fs);
            end

            T = repmat(tindx,1,length(wint))+repmat(wint,length(tindx),1)+1;
            T(T<1) = 1;
            T(T>length(xin))=length(xin);

%                     for k = 1:length(xin)
                Xchop = xin(T);

                hasnans = any(xisnan(T));
                segment.discarded = hasnans;
                if all(hasnans) && all(islogical(makeplot)) && makeplot>=0
                    fprintf('\nAll segments contain NaN values. Discarding these data')
                    if nargout > 1
                         varargout = {Xsh,Xwin,T,wint,segment};
                    end
                    return
                elseif any(hasnans) && ~any(ishandle(makeplot)) && makeplot>=0
                    fprintf('\n%i (%0.2f %%) Segments with NaN values have been excluded',sum(hasnans),100*mean(hasnans))
                    T = T(:,~hasnans);
                    segment.wint = segment.wint(~hasnans);
                    Xchop = xin(T);
                end

%                     end
            if length(xin)==1  %%% Streaming is only implemented for single-channel input
                snip = xin(T(end)+1:numel(xin));
                if ~isempty(snip) && ~any(isnan(snip))                        
                    me(1).write_buffer(snip);
                end
            end
        else
            Xchop = xin;
            T=0;
            wint=[];
        end
%                 Xchop(end+1:me(1).fftN,:) = 0;
        del = Inf;
%                 tol =1; % Stop when the average shift is less than 1 sample
        tol =me(1).sampling_rate/me(1).lowpass(1);
        k = 0;
        olddt2 = 0;
        olddt =0;
%                  Xwin = Xchop.*repmat(kaiser(me(1).bufferN,5),1,size(Xchop,2));
%                 for k = 1:length(Xchop)
        Xwin = Xchop.*repmat(me(1).win,1,size(Xchop,2));
%                    Xwin = Xchop.*repmat(me(1).win,1,size(Xchop,2)); 
        Xwin(end+1:me(1).fftN,:)=0;
        Xsh = Xwin;
        Xfilt=Xsh;
%                 end
        if   double( makeplot(1))>=0                     
            fprintf('\nComponent %3i Iter %3i',compno,0)
        end
        color_cycle = 10;
        if all(ishandle(makeplot))
            set(makeplot(1:end-1),'ydata',zeros(me(1).bufferN,1));
        end
%                 std_moment = @(x)mean(nanmean(x.^me(1).order)./(nanmean(x.^2)).^(me(1).order/2));
        if isempty(me(1).sampweight)
            smpw = 1;
        else
            smpw = me(1).sampweight;
        end
%                 std_moment = @(x)mean(cumulant(x,me(1).order,1,smpw')./(nanmean(x.^2).*nanmean(smpw.^2)).^(me(1).order/2));
%         std_moment = @(x)mean(cumulant(x,me(1).order,1)./(nanmean(x.^2).*nanmean(smpw.^2)).^(me(1).order/2));
        std_moment = @(x)mean(cumulant(x,me(1).order,1,true));
         switch me(1).order
            case 3
                moment_type = 'skewness';
            case 4
                moment_type = 'kurtosis';
            otherwise
                moment_type = 'standardized cumulant';
         end
        apply_window = false;
        use_shifted=false;
%                 initialize = true;
        if me(1).use_partial_delay_method
            updatebsp = me(1).do_bsp_update;
            me(1).do_bsp_update = false; %Avoid updating the bispectrum estimate twice unnecessarily here.
            me(1).get_input(Xsh,apply_window,use_shifted,initialize);
            me(1).do_bsp_update = updatebsp;
            [Gpart,sgn] = me(1).partial_delay_filt(Xsh,false,true); %
            if ~isequal(sgn,1)
                Xsh = Xsh.*sgn;
                Xwin = Xwin.*sgn;
            end
            Gpart(isnan(Gpart))=0;
            if me(1).do_filter_update
                me(1).G = mean(Gpart,2);
            end
%                     me(1).apply_filter(Xsh,false,true);
%                     newdt = me(1).delay;
        end
        while del >tol && k < maxiter                    
            if all(ishandle(makeplot))
                set(makeplot(1),'cdata',Xsh');

                set(makeplot(7),'string',sprintf('Component %3i, Iter. %3i, Mean shift = %2.2fs, %s=%2.2f',compno,k,del/me(1).sampling_rate,moment_type,std_moment(Xfilt)));
                %                         set(makeplot(mod(k,5)+2),'ydata',ifftshift(abs(me(1).filterfft)));
%                         for pli = 1:length(makeplot)
%                             set(plh(pli),'ydata',get(makeplot(pli),'ydata')+std(me(1).feature)*.5);
%                          end
                set(makeplot(mod(k,5)+2),'ydata',me(1).feature,'Color',hsv2rgb([mod(k,color_cycle)/color_cycle 1 .8]));
                drawnow
            elseif double(makeplot) > 0
                figure,
                subplot(2,1,1)
                makeplot = imagesc(fftshift(me(1).sampt)/me(1).sampling_rate,[],Xsh');
                makeplot(7) = title(sprintf('Component %03i, Iter. %3i, Mean shift = %2.2fs, %s=%2.2f',compno,k,del/me(1).sampling_rate,moment_type,std_moment(Xfilt)));
                subplot(2,1,2)
%                          makeplot(2:6) = plot(ifftshift(me(1).freqs),ifftshift(abs(me(1).filterfft))*ones(1,5));
                 plh = plot(fftshift(me(1).sampt)./me(1).sampling_rate,me(1).feature*ones(1,5));
%                          cmap = hsv(length(plh));
                 for pli = 1:length(plh)

                     set(plh(pli),'Color',hsv2rgb([mod(pli,color_cycle)/color_cycle 1 .8]));
                 end
                 makeplot(2:6)=plh;
%                         xlim([0 me(1).lowpass])
                drawnow
            end

            k=k+1;
            if  double( makeplot(1))>=0
                fprintf('\b\b\b%03i',compno,k)
            end

            if k >2
              olddt = me(1).delay;
            end

            

           if me(1).use_partial_delay_method
               [Xfilt,FXsh,sgn] = me(1).apply_filter(Xsh,false,true);
%                        Xsh = real(ifftshift(ifft(FXsh),1));
               Xsh = real(ifft(FXsh));
               newdt = me(1).delay;
               delt = me(1).radw(me(1).keepfreqs{1})*newdt;
               Gpart = (sgn.*exp(-1i.*delt)).*Gpart; 
               Gpart(isnan(Gpart))=0;
               G = mean(Gpart,2);
               if me(1).do_filter_update
                   me(1).G = G;
               end
               if me(1).do_wave_update
                me(1).feature = nanmean(Xsh,2);
               end
              
               if me(1).adjust_lag && me(1).do_filter_update
%                            ffun = ifftshift(real(ifft(me(1).filterftlag)));                   
                   ffun = ifftshift(real(ifft(me(1).filterftlag.*abs(me(1).waveftlag+eps))),1);   
                   mph = sum(exp(-1i*2*pi*me(1).sampt(:)./me(1).fftN).*abs(ffun).^2)./sum(abs(ffun).^2);                   
                   mph = mph./(abs(mph)+eps);
                   if ~isnan(mph)
                     me(1).lag = mph; % Circular shift to keep filter energy centered on the window
                   end

               end

           else
               %%% This approach recomputes all statistics at every
               %%% iteration, which is unnecessary.
                me(1).get_input(Xsh,apply_window,use_shifted,initialize);
                [Xfilt,FXsh] = me(1).apply_filter(Xsh,false,true);
               newdt = me(1).delay;
%                         Xsh = real(ifftshift(ifft(FXsh),1));
                Xsh = real(ifft(FXsh));
           end
            % checks two and one step back to reduce getting
            % trapped at points of cyclical stability.
%                     del = sqrt(mean((olddt2-newdt).^2));%min(sqrt(mean((olddt-newdt).^2)),sqrt(mean((olddt2-newdt).^2)));
            del = std(olddt2-newdt);%min(sqrt(mean((olddt-newdt).^2)),sqrt(mean((olddt2-newdt).^2)));
            olddt2 = olddt;


%                     if ~isempty(me(1).regval)
%                     %%% Now update the regressor weights, if there are any
%                     %%% regressors
%                     
%                         cumXfilt = cumulant(Xfilt,me(1).order);
%                         nrm = @(x)x./sqrt(sum(abs(x).^2));
% %                         me(1).regweight = nrm(me(1).regval.value\cumXfilt');
% %                         me(1).sampweight=nrm(me(1).regval.value*me(1).regweight);
%                     end
          
           
  
           
        end
        
        %%% With highly periodic signals in high noise,
        %%% the algorithm seems occasionally to converge with inverted sign (for odd orders).
        %%% To avoid this, make sure that partial polycoherence is
        %%% positively correlated with sample polycoherence.
        [Xfilt,FXsh] = me(1).apply_filter(Xwin,false,true);
        if mod(me(1).order,2)~=0 && skewness(Xfilt(:))<0 %me(1).pbcorr < 0 
           fprintf('\nInverting sign...')
           me(1).waveftlag = -me(1).waveftlag;
           me(1).G= -me(1).G;
           [~,FXsh] = me(1).apply_filter(Xwin,false,true);
        end
        
        %%% Also need to make sure that the output of the filter applied to
        %%% the feature waveform is centered with respect to the maximum!
     
         Xsh = real(ifft(FXsh));
         me(1).feature = nanmean(Xsh,2);
        [~,mxi] = max(real(ifft(me(1).filterftlag.*me(1).waveftlag+eps)));
        if mxi~=1 && ~isnan(mxi) && me(1).do_filter_update
            me(1).filterfun = circshift(me(1).filterfun,ceil(-me(1).sampt(mxi)/2));
            me(1).feature= circshift(me(1).feature,floor(-me(1).sampt(mxi)/2));
        end
        
       if all(ishandle(makeplot))
            set(makeplot(1),'cdata',Xsh');
            set(makeplot(6),'ydata',me(1).feature,'Color','k','linewidth',2);
            drawnow
       end
       
       
        if me(1).power_iterate && me(1).do_filter_update && me(1).do_wave_update
                 for rep = 1:50
                    me(1).power_iteration(Xwin);
                 end
        end
        %%% Set the delays to the correct value for the original
        %%% data set;

        [Xfilt,~] = me(1).apply_filter(Xwin,apply_window);
         me(1).running_mean = nanmean(Xfilt(:));
          me(1).running_ssq = nanmean(Xfilt(:).^2);
          me(1).running_var = me(1).running_ssq-me(1).running_mean.^2;
%                   Xfilt = (Xfilt-me(1).running_mean)./sqrt(me(1).running_var);
%                  Xfilt = (Xfilt-nanmean(Xfilt))./nanstd(Xfilt);

        if ~isempty(segment.wint)
            segment.wintadj = segment.wint+me(1).delay/segment.fs;
        end
        %%% Update CDF buffer
            xfilt = (me(1).xfilt(xin)-me(1).running_mean)./sqrt(me(1).running_var);
        if size(xin,1)>me(1).buffersize
            Tcdf = chopper([0 me(1).buffersize*me(1).CDFupsample-1],segment.wint(1:me(1).CDFupsample:end),segment.fs);
%                     Tcdf(Tcdf<1)=1;Tcdf(Tcdf>length(xin))=length(xin);
            Tcdf = mod(Tcdf-1,length(xfilt))+1;
            Xfilt = sort(xfilt(Tcdf));
        else
            nw = floor(numel(xfilt)./(size(xfilt,1)*me(1).CDFupsample));

            Xfilt = reshape(xfilt(1:nw*(size(xfilt,1)*me(1).CDFupsample)),(size(xfilt,1)*me(1).CDFupsample),nw);
            Xfilt = sort(Xfilt);
        end

       if me(1).use_adaptive_threshold
           for k = 1:me(1).threshold_order
               me(1).CDFbuffer(:,k) = nanmean(Xfilt.^k,2);
           end
       end    

        if isempty(segment.wint)
%                     T = mod(repmat((0:size(Xsh,1)-1)',1,size(Xsh,2))+repmat(me(1).delay,size(Xsh,1),1),size(Xsh,1))+1;
%                     T = T+ones(size(T,1),1)*(0:size(T,2)-1)*size(T,1);
        else
            T = T+ repmat(me(1).delay,size(T,1),1);
            T(T<1)=1;
            T(T>length(xin))=length(xin);

            segment.wintadj = segment.wint+me(1).delay/segment.fs;
        end
    else
        me(1).write_buffer(xin);
        T=[];
    end    
    me(1).segment = segment;
  
    if length(me)>1
       xrec = me(1).reconstruct(xin);
       if nargout > 0
           [Xsh2,Xwin2,T2] = me(2:end).get_block(xin-xrec,maxiter,makeplot,segment,compno+1,initialize);
           Xsh = cat(3,Xsh,Xsh2);
           Xwin = cat(3,Xwin,Xwin2);
           T = cat(3,T,T2);
       else
           if  length(me(2).B)==length(me(1).B)
                me(2).Imats = me(1).Imats;
                me(2).Iconjmats = me(1).Iconjmats;                    
           end
            me(2:end).get_block(xin-xrec,maxiter,makeplot,segment,compno+1,initialize);
       end
    end
    if nargout > 0
        varargout = {Xsh,Xwin,T,wint,segment,makeplot};
    end
end