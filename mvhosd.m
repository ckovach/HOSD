
classdef mvhosd < hosobject
    
    properties
     
      Gpart
      Xwin
      maxplotn =4;
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
            
            std_moment = @(x)mean(cumulant(x,me(1).order,1)./(nanmean(x.^2).*nanmean(smpw.^2)).^(me(1).order/2));
             switch me(1).order
                case 3
                    moment_type = 'skewness';
                case 4
                    moment_type = 'kurtosis';
                otherwise
                    moment_type = 'standardized cumulant';
             end
             nsig = size(Xsh,3);
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

               [Xfilt,FXsh,sgn] = me(:,1).apply_mvfilter(Xsh,false,true);
               Xsh = real(ifft(FXsh));
               
               newdt = me(1).delay;
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
                [~,mxi] = max(sgn.*L);
                flcorrection = me(1).sampt(mxi);
                delt2 = me(1).radw*flcorrection; 
                Gpart = Gpart.*exp(1i*permute(delt2,[1 3 2]));
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
                FXfilt(me(1).keepfreqs{1},:) = sum(FXwin(me(1).keepfreqs{1},:,:).*repmat(me(1).filterfft(me(1).keepfreqs{1},:,:),1,size(X,2)),3);
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
                filts = squeeze(me(1).filterfun);
                Xin = X;
                Xin(end+me(1).fftN,:) = 0;
                Xfilt = 0;
                if nargout > 3
                    mvXfilt = zeros(size(Xin,1),size(filts,2));
                end
                for k = 1:size(filts,2)
                    xf = filter(filts(:,k),1,Xin(:,k));
                    Xfilt = Xfilt+xf;
                    if nargout > 3
                        mvXfilt(:,k) = xf;
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
           if size(in,chdim) ~= size(me(1).feature,3) && size(in,chdim-1) == size(me(1).feature,3)
               warning('MVHOS expected dimension %i for channels and %i for features, but size suggests they are reversed.\nThese will be exchanged now. In the future make sure the dimensions are correctly ordered,\nas this would have been missed if the number of features and channels happened to coincide.',chdim,chdim-1)
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
                   mvout = cat(sum(size(in)>1)+1,mvout,mvxf);
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
           out = me(1).reconstruct(in,thresh,apply_window,varargin{:}); 
           if length(me)>1
%                if size(in,2) == 1
%                    applydim = 2;
%                else 
%                    applydim = max(find(size(in)>1))+1;
%                end
               out =  cat(sum(size(in)>1)+1,out,me(2:end).xrec(in-out(:,1),thresh,apply_window,varargin{:}));
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

        
    end

end



