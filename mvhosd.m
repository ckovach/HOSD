
classdef mvhosd < hosobject
    
    properties
     
        Gpart
        Xwin
        maxplotn =4;
      
        %%% If annealing_start is greater than 0, then Gaussian white noise
        %%% will be added to the input during iterated realignment to improve convergence
        %%% and  decremented with iteration according to annealing_schedule.
        annealing_start=0;%Starting noise amplitude used for annealing, in units of input s.d
        annealing_schedule = @(k,maxk)((maxk-k)/maxk); %How to scale annealing noise as a function of iteration number (1st arg.) and maximum iterations (2nd arg)    

    end
    
    methods

        function me = mvhosd(varargin)

            
            me = me@hosobject(varargin{:});
            
        end
        
        function [Xsh,Xwin,makeplot] = align(me,Xwin,Gpart,maxiter,makeplot,compno)
             % Fit a block of data all at once
            % Process input if length is >= buffer size, else add to buffer.
            if nargin < 6
                compno = 1;
            end
            if nargin < 5
                makeplot = true;
            end
            if nargin < 4 || isempty(maxiter)
                maxiter = 25;
            end
%             if ~iscell(xin)
%                 xin = {xin};
%             end
%             nxin = numel(xin);
%             xisnan = isnan(xin);       
            
            if size(Xwin,3)>1
               nsig = size(Xwin,3);
               sigdim = 3;
            else
                nsig = size(Xwin,2);
                sigdim = 2;
            end
            
            me(1).bufferPos = 0; % Discard the buffer
           
%             Gpart = me(1).Gpart;
            
            Gpart0 = Gpart;
            
      %      Xwin(end+1:me(1).fftN,:) = 0;
%             Xwin = me(1).Xwin;
            Xsh = Xwin;
             Xfilt=nanmean(Xsh,3);
           
                 
            del = Inf;
            tol =me(1).sampling_rate/me(1).lowpass(1);
            iter = 0;
            olddt2 = 0;
            olddt =0;
          
            fprintf('\nComponent %3i Iter %3i',compno,0)
            color_cycle = 10;
            if all(ishandle(makeplot))
                set(makeplot(1:end-1,:),'ydata',zeros(me(1).bufferN,1));
            end
            
            if isempty(me(1).sampweight)
                smpw = 1;
            else
                smpw = me(1).sampweight;
            end
            
%             std_moment = @(x)nanmean(cumulant(x,me(1).order,1)./(nanmean(x.^2).*nanmean(smpw.^2)).^(me(1).order/2));
%             std_moment = @(x)nanmean(cumulant(x,me(1).order,1)./(nanmean(x.^2)).^(me(1).order/2));
            std_moment = @(x)nanmean(cumulant(x,me(1).order,1,true));
            switch me(1).order
                case 3
                    moment_type = 'skewness';
                case 4
                    moment_type = 'kurtosis';
                otherwise
                    moment_type = 'standardized cumulant';
            end
            nsig = size(Xsh,3);
            if me(1).annealing_start>0
                noise = randn(size(Xsh)).*nanstd(reshape(Xsh,[size(Xsh,1)*size(Xsh,2),1,size(Xsh,3)]));
            else
                noise = 0;
            end
            
            while del >tol && iter < maxiter                  
                [~,plotcompi] = sort(sum(abs(me(1).wavefft).^2.*abs(me(1).filterfft).^2),'descend');
                
                try
                     if all(ishandle(makeplot))
                       for k = 1:size(makeplot,2)
                            kplot = plotcompi(k);
                            set(makeplot(1,k),'cdata',Xsh(:,:,kplot)');

                            set(makeplot(7,k),'string',sprintf('Ch.%3i, Comp.%3i, Iter.%3i\nMean shift =%2.2fs, %s=%2.2f',kplot,compno,iter,del/me(1).sampling_rate,moment_type,std_moment(Xfilt)));
                            set(makeplot(mod(iter,5)+2,k),'ydata',me(1).feature(:,kplot),'Color',hsv2rgb([mod(iter,color_cycle)/color_cycle 1 .8]));
                            axis tight
                            ylim(minmax(me(1).feature(:)));
                        end
                        drawnow
                    elseif islogical(makeplot) && makeplot
                        figure,
                        clear makeplot
                        nplot = min(nsig,me(1).maxplotn);
                        for k = 1:nplot
                            kplot = plotcompi(k);
                            subplot(2,max(nplot,me(1).maxplotn),k )
                            makeplot(1,k) = imagesc(fftshift(me(1).sampt)/me(1).sampling_rate,[],Xsh(:,:,kplot)');
                            makeplot(7,k) = title(sprintf('Ch.%3i, Comp.%3i, Iter.%3i\nMean shift =%2.2fs, %s=',kplot,compno,iter,del/me(1).sampling_rate,moment_type ));
                            subplot(2,nplot,k+ nplot)
                           plh = plot(fftshift(me(1).sampt)./me(1).sampling_rate,me(1).feature(:,kplot)*ones(1,5));
                             for pli = 1:length(plh)

                                 set(plh(pli),'Color',hsv2rgb([mod(pli,color_cycle)/color_cycle 1 .8]));
                             end
                             makeplot(2:6,k)=plh;
                        end
    %                         xlim([0 me(1).lowpass])

                     end
                    drawnow
                catch
                        plh = false;
                end
                iter=iter+1;
                fprintf('\b\b\b%03i',compno,iter)

                if iter >2
                  olddt = me(1).delay;
                end
                if me(1).annealing_start>0
                    Xsh0 = Xsh;
                     
%                     noise = randn(size(Xsh))*nanstd(Xsh(:));
             
                    Xsh = Xsh + noise*me(1).annealing_start*me(1).annealing_schedule(iter,maxiter);
                    [Xfilt,~,sgn] = me(:,1).apply_mvfilter(Xsh,false,true);
                    newdt = me(1).delay;
                    FXsh = fft(Xsh0).*conj(fft(me(1).sampt'==newdt)); 
                else
                    [Xfilt,FXsh,sgn] = me(:,1).apply_mvfilter(Xsh,false,true);
                    newdt = me(1).delay;
                end
                Xsh = real(ifft(FXsh));
               
%                delt = me(1).radw(me(1).keepfreqs{1})*newdt;
               delt = me(1).radw*newdt;

               %                    
%                if isscalar(smpw)
%                     [~,mxi] = max(Xfilt.^me(1).order);
%                 else
%                     [~,mxi] = max(Xfilt.^me(1).order.*repmat(smpw',size(Xfilt,1),1));
%                 end
%                 if mod(me(1).order,2)==0 
%                     sgn = sign(Xfilt(mxi + (0:size(Xfilt,2)-1)*size(Xfilt,1)).*smpw'); 
%                 end
%                 
%                 samptc=(me(1).sampt);
%                 dt = samptc(mxi);
%                 delt = me(1).radw*dt;
%                 newdt = dt;
               %  deltg = delt(me(1).keepfreqs{1},:);
                Gpart = (sgn.*exp(-1i.*delt)).*Gpart; 
                G = nanmean(Gpart,2);
           
%                Xsh = real(ifftshift(ifft(exp(1i*delt).*FXwin),1));
%                Xsh = real(ifft(exp(1i*delt).*FXwin));
                features = nanmean(Xsh,2);
           
                %%% G should be matched to features so as to produce a
                %%% peak at zero lag. Might as well enforce this
                %%% explicitly as aligment arrors might accrue.
                
                L = ifft(G.*fft(features));
                [~,mxi] = max(real(L).^(2-mod(me(1).order,2)));
                
                if mod(me(1).order,2)==0
                    sgc=sign(real(L(mxi+permute((0:size(L,3)-1)*size(L,1),[1 3 2]))));
                else
                    sgc=1;
                end
                flcorrection = me(1).sampt(mxi);
                delt2 = me(1).radw*flcorrection; 
                Gpart = Gpart.*exp(1i*permute(delt2,[1 3 2])).*sgc;
                G = mean(Gpart,2);
                
                if size(G,1) == me(1).fftN
                    G = G(me(1).keepfreqs{1},:,:);
                end
                me(1).G = G;
                me(1).waveftlag= fft(features);
 
                if me(1).adjust_lag
                   ffun = ifftshift(sum(real(ifft(me(1).filterftlag.*abs(me(1).waveftlag+eps))),3),1);                   
                   mph = sum(exp(-1i*2*pi*me(1).sampt(:)./me(1).fftN).*sum(abs(ffun).^2,2))./sum(sum(sum(abs(ffun).^2,2),3));                   
                   mph = mph./(abs(mph)+eps);
                   me(1).lag = mph; % Circularshift to keep filter energy centered on the window
                end
                Gpart = Gpart;
                Xwin = Xsh;
              
            
                % checks two and one step back to reduce getting
                % trapped at points of cyclical stability.
                del = std(olddt2-newdt);%min(sqrt(mean((olddt-newdt).^2)),sqrt(mean((olddt2-newdt).^2)));
                olddt2 = olddt;


            end

            %%% Set the delays to the correct value for the original
            %%% data set;
             [~,~] = me(1).apply_mvfilter(Xwin,false,true);
%             T = T+ repmat(me(1).delay,size(T,1),1);
%             T(T<1)=1;
%             T(T>length(xin))=length(xin);
            if size(me,2)>1
               Xrec = me(:,1).xrec(Xwin);
%                if nargout > 0
%                    [Xsh2,Xwin2,T2] = me(2:end).get_block(xin-xrec,maxiter,makeplot,compno+1);
%                    me(2).Xwin = Xwin-Xrec;
                   delt = me(1).radw*me(1).delay; 
%                    me(2).Gpart = (Gpart - me(1).filterfft).*exp(1i.*delt); % The partial delay filter has to retain the correct alignment to Xwin.
%                    [Xsh2,Xwin2] = me(2:end).get_block(xin-xrec,maxiter,makeplot,compno+1);
                   [Xsh2,Xwin2] = me(2:end).align(Xwin-Xrec,(Gpart - me(1).filterfft).*exp(1i.*delt),maxiter,makeplot,compno+1);
                   Xsh = cat(3,Xsh,Xsh2);
                   Xwin = cat(3,Xwin,Xwin2);
%                    T = cat(3,T,T2);
%                else
%                     me(2:end).get_block(xin-xrec,maxiter,makeplot,compno+1);
%                end
            end
        end
        
        %%%%%%%%%
        function [Xfilt,FXshift,sgn,mvXfilt] = apply_mvfilter(me,X,apply_window,return_shifted,varargin)
            if nargin<3 || isempty(apply_window)
                apply_window = true;
            end
            if nargin < 4 || isempty(return_shifted)
               return_shifted = true; 
            end           
            FXshift = [];
            sgn = 1;
            if isempty(me(1).sampweight)
               smpw =1 ;
            else
                smpw = me(1).sampweight;
            end
            if size(X,1) == me(1).bufferN
                if apply_window
                    win = me(1).win;
                else
                    win =ones(size(X,1),1);
                end
               % Xwin = fftshift(repmat(win,1,size(X,2),size(X,3)).*X,1);
                Xwin = repmat(win,1,size(X,2),size(X,3)).*X; %#ok<*PROPLC>
                Xwin(end+1:me(1).fftN,:,:) = 0;
                FXwin = fft(Xwin);
%                 FXwin = FXwin(me(1).keepfreqs{1},:,:);
%                filts = cat(3,me(:,1).filterfft);
                
%                 FXwin = fft(X)';
                FXfilt =zeros(size(Xwin,1),size(Xwin,2));
                FXfilt(me(1).keepfreqs{1},:) = nansum(FXwin(me(1).keepfreqs{1},:,:).*repmat(me(1).filterfft(me(1).keepfreqs{1},:,:),1,size(X,2)),3);
                Xfilt = real(ifft(FXfilt));
                
%                 Xfilt = sum(real(ifft(FXwin.*repmat(me(1).filterfft(me(1).keepfreqs{1}),1,size(X,2)))),3);   
                if nargout > 3
                    mvXfilt = real(ifft(FXwin.*repmat(me(1).filterfft,1,size(X,2))));   
                end                    
                if isscalar(smpw)
                    [~,mxi] = max(Xfilt.^me(1).order);
                else
                    [~,mxi] = max(Xfilt.^me(1).order.*repmat(smpw',size(Xfilt,1),1));
                end
                if mod(me(1).order,2)==0 
                   sgn = sign(Xfilt(mxi + (0:size(Xfilt,2)-1)*size(Xfilt,1)).*smpw'); 
                end
                if nargout >1 && return_shifted
                 
    %                 FX = fft(X);
                    samptc=(me(1).sampt);
                     dt = samptc(mxi);
                     me(1).delay = dt;             
                    
                    delt = me(1).radw*dt;
                    FXshift = exp(1i*delt).*FXwin;
       
                else
                     FXshift = FXwin;
                end
                FXshift = FXshift*diag(sgn);
            else
%                  Xin = X(:);
%                 filts = squeeze(me(1).filterfun);
                filts = me(1).filterfun;
                if size(X,3) ~=size(filts,3)
                    X = permute(X,[1 3 2]);
                end
                Xin = X;
                Xin(end+me(1).fftN,:) = 0;
                Xfilt = 0;
                if nargout > 3
                    mvXfilt = zeros(size(Xin,1),1,size(filts,3));
                end
                for k = 1:size(filts,3)
                    xf = filter(filts(:,:,k),1,Xin(:,:,k));
                    Xfilt = Xfilt+xf;
                    if nargout > 3
                        mvXfilt(:,:,k) = xf;
                    end
                end
                Xfilt = Xfilt(ceil(me(1).fftN/2)+1:end-floor(me(1).fftN/2),:);
                if nargout > 3
                    mvXfilt = mvXfilt(ceil(me(1).fftN/2)+1:end-floor(me(1).fftN/2),:,:);
                end
            end
               
        end
 %%%%%%%
        function [out,mvout] = xfilt(me,in,apply_window)
           if nargin < 2
               in = me(1).dat;
           end
           if nargin < 3 || isempty(apply_window)
              apply_window = false; 
           end
           chdim = find(size(me(1).feature)>1,1,'last');
           if ~isvector(in) && size(in,chdim) ~= size(me(1).feature,3) && size(in,chdim-1) == size(me(1).feature,3)
%                warning('MVHOS expected dimension %i for channels and %i for features, but size suggests they are reversed.\nThese will be exchanged now. In the future make sure the dimensions are correctly ordered,\nas this would have been missed if the number of features and channels happened to coincide.',chdim,chdim-1)
               in = permute(in, [1:chdim-2 chdim chdim-1]);
           end
           if nargout > 1
               [out,~,~,mvout] = me(1).apply_mvfilter(in,apply_window,false);  
           else
               [out,~] = me(1).apply_mvfilter(in,apply_window,false);
           end
%             if size(in,1)==me(1).buffersize %% Make sure the output is consistent if the input happens to be of buffersize length
%                 out = ifftshift(out,1);
%             end
            
           if size(me,2)>1
               if nargout == 1
                   xf = me(:,2:end).xfilt(in-me(1).xrec(in,[],apply_window),apply_window);
               else                   
                   [xf,mvxf] = me(:,2:end).xfilt(in-me(1).xrec(in,[],apply_window),apply_window);
                   mvout = cat(find([size(in),1]==1,1),mvout,mvxf);
               end
               out = cat(sum(size(in)>1),out,xf);
           end
        end
        %%%%%%%
       
        function [Xrec,Xfilt] = reconstruct(me,X,threshold,apply_window)
            
            if nargin < 2
                X = [me(:,k).inputbuffer];
            end
            if nargin < 3 || isempty(threshold)
                threshold = me(1).thresh;
            end
            if nargin < 4 
                apply_window = [];
            end
            xisnan = isnan(X);
            Xfilt = me(1).apply_mvfilter(X,apply_window);
             if size(X,1) == me(1).bufferN
                 Xfilt = ifftshift(Xfilt,1);             
             end
            if size(X,3)>1
                applydim = 3;
            else
                applydim = 2;
            end

            Xthr=me(1).filter_threshold(Xfilt,threshold);
            
        %     wf = fftshift(me(1).waveform,1);
            wf = me(1).waveform;
            wf(end+1:size(Xthr,1),:) = 0;
            wf = circshift(wf,-floor(me(1).bufferN/2));
            Xthr(end+1:length(wf),:)=0;
            FXthresh =fft(Xthr);
            featfft = fft(wf);          
            if applydim == 3
               featfft=repmat(featfft,1,size(X,2));
            end
            Xrec = real(ifft(FXthresh.*featfft));
            Xrec(size(X,1)+1:length(wf),:) = [];
            %X(xisnan)=0;
            Xrec(xisnan)=0;

            use_filtered_lmse = true;
            if use_filtered_lmse
                %%% Apply the filter to the reconstructed data for LMSE fitting
                %%% so that the frequencies are appropriately weighted.
                Xrecfilt = me(1).xfilt(Xrec,apply_window);
                if size(X,1) == me(1).bufferN
                     Xrecfilt = ifftshift(Xrecfilt,1);
                end
                
                Xrecfilt(isnan(Xfilt))=0;
                Xfilt(isnan(Xfilt)) = 0;
                
                beta = Xrecfilt(:)'*Xfilt(:)./sum(Xrecfilt(:).^2);
                Xrec = beta*Xrec; 
            else
                a= sum(abs(Xrec(:)).^2); %#ok<*UNRCH>
                if a > 0
                 Xrec = Xrec*(X(:)'*Xrec(:))./a; % Scale to minimize total mse.
                end
            end

            Xrec(xisnan) = nan;
%             if nargin < 2
%                 me(1).reconbuffer = Xrec;
%             end
        end
        
       
        %%%%%%%
         function out = xthresh(me,in,threshold,apply_window)
           if nargin < 2
               in = me(1).dat;
           end
            if nargin < 3 || isempty(threshold)
                threshold = me(1).thresh;
            end
            if nargin < 4 || isempty(apply_window)
                apply_window = false;
            end
            Xfilt = me(1).xfilt(in);
            out=me(1).filter_threshold(Xfilt,threshold);
            if size(me,2)>1
               out = cat(sum(size(in)>1),out,me(2:end).xthresh(in-me(1).xrec(in,threshold,apply_window),threshold,apply_window));
           end
         end
           %%%
        function out = xrec(me,in,thresh,apply_window,varargin)
           if nargin < 2
               in = me.dat;
           end
           if nargin < 3
               thresh = [];
           end
           if nargin < 4
               apply_window = false;
           end
           if size(in,3) == 1
               in = permute(in,[1 3 2]);
           end
           
           out = me(1).reconstruct(in,thresh,apply_window,varargin{:}); 
           if length(me)>1
%                if size(in,2) == 1
%                    applydim = 2;
%                else 
%                    applydim = max(find(size(in)>1))+1;
%                end
               out =  cat(find([size(in),1]==1,1),out,me(2:end).xrec(in-out,thresh,apply_window,varargin{:}));
           end
        end
         %%%%
          function out = ximp(me,in,apply_window)
           if nargin < 2
               in = me(1).dat;
           end
           if nargin < 3 || isempty(apply_window)
              apply_window = false; 
           end
            
            Xthr=me(1).xthresh(in);
            out = Xthr>0;
            if size(me,2)>1
               out = cat(sum(size(in)>1),out,me(:,2:end).ximp(in-me(:,1).xrec(in,[],apply_window),apply_window));
            end
          end

          function [A,B,makeplot] = get_block(me,in,maxiter,makeplot,segment,compno,initialize)
            
              if nargin < 6 || isempty(compno)
                  compno = 1;
              end
              if nargin < 5 || isempty(segment)
                  segment = [];
              end
              if nargin < 7 || isempty(initialize)
                  initialize = true;
              end
              if nargin < 4 || isempty(makeplot)
                  makeplot = true;
              end
              if nargin < 3 || isempty(maxiter)
                  maxiter = 25;
              end
              if size(in,2) > 1 && size(in,3) ==1
                   in = permute(in,[1 3 2]);
              end
               
              if all(isnan(in))
                  A = [];
                  B = [];
                  return
              end
                   
                nfp = 0;
                nchan = size(in,3);
                for k = 1:nchan
                    nfp=fprintf([repmat('\b',1,nfp),'\nComp. %i,estimating HOS for chan. %i'],compno,k)-nfp;
                    X(:,:,k) = me(1).chop_input(in(:,k));
                    Gpart(:,:,k) = me(1).partial_delay_filt(X(:,:,k),true,true); 
                    if k==1
                        X(:,:,nchan)=0;
                        Gpart(:,:,nchan)=0;
                    end
                end


                me(1).use_adaptive_threshold=false;
                
                me(1).filterfft = nanmean(Gpart,2);
                me(1).feature = nanmean(X,2);
                [A,B,makeplot] = me(1).align(X,Gpart,maxiter,makeplot,compno);
                if size(A,4) ==1
                    A = permute(A,[1 2 4 3]);
                end
%                  [~,xfilt] = me(1).xfilt(permute(zresid,[1 3 2]));
               
%                 if nk > 1 && skf(nk-1)< params.skewness_threshold && skf(nk)< params.skewness_threshold 
%                     fprintf('Skewness under %0.3f for the last 2 components... stopping at %i.',params.skewness_threshold,nk)
%                     return
%                 else
                if length(me)>1
                    xrec = me(1).xrec(in);
                     [A2,B2] = me(2:end).get_block(in-xrec,maxiter,makeplot,segment,compno+1,initialize);
                    A = cat(3,A,A2);
                    B = cat(3,B,B2);
                end
              
          end
        
    end

end



