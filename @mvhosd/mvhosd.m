
classdef mvhosd < hosobject
    
    % Multivariate HOSD estimation using a auto-HOS and (optionally) static ICA using a
    % power-iteration technique (PICA). MVHOSD interleaves the estimation of a multivariate filter
    % (i.e. a spatio-temporal filter) and re-weightng of the spatial
    % dimension with static ICA estimates. Both steps maximize the cumulant
    % of a given order. 
    %
    % USAGE:
    %
    % To define an MVHOSD object, first create a univariate HOSOBJECT with
    % the desired order and component number, then instantiate the MVHOSD
    % object using the constructor:     
    %       hos = hosobject(...);
    %       mvhos = mvhosd(hos);
    %Then run the decomposition on data, X:
    %       mvhos.get_block(X);
    % where X is a time x channel OR time x trial x channel matrix. 
    %
    % OUTPUT:
    % 
    % After estimation, to obtain a time x component x channel matrix of features, F, use
    %       F = [mvhos.feature];
    % and the a matrix representing a time x component x channel
    % spatiotemporal unmixing/deconvolution filter is obtained with
    %       H = [mvhos.filterfun];
    %
    % To obtain the filtered (deconvolved/unmixed) signal, with serial deflation of each
    % component:
    %       [xfilt,XF] = mvhos.xfilt(X);
    % where XF is separated by channel and xfilt = sum(XF,2).
    %
    % To obtain the filtered and thresholded filter use:
    %       xthresh = mvhos.xthresh(X);
    % To obtain the time x channel x component reconstruction:
    %       xrec = mvhos.xrec(X);
    %
    % See also HOSOBJECT and PICA
    
    %C. Kovach 2019-2022
    
    properties
     
        Gpart
        Xwin
        maxplotn =4;
      
        %%% If annealing_start is greater than 0, then Gaussian white noise
        %%% will be added to the input during iterated realignment to improve convergence
        %%% and  decremented with iteration according to annealing_schedule.
        annealing_start=0;%Starting noise amplitude used for annealing, in units of input s.d
        annealing_schedule = @(k,maxk)((maxk-k)/maxk); %How to scale annealing noise as a function of iteration number (1st arg.) and maximum iterations (2nd arg)    
    
        %%% Interleave filter estimation with static ICA on the filter
        %%% output with the hope of improving spatial separation.
        do_static_ica = 10; %Updates every kth iteration (0 = none)
       
        subspace_dim = 0; %Dimensionality of the feature. The response in channel i, freq. f is modeled as Hi(f) = Ai*B'(f).
                          %Subspace_dim is the rank of Ai and B. If the rank is 1, each channel is assumed to contain a 
                          %scaled copy of the same waveform. Subspace_dim=0 is full rank (no dimensionality reduction).
        projection = 1;   %Current static projection if subspace_dim > 0;
      end
    
    methods

        function me = mvhosd(varargin)

            
            me = me@hosobject(varargin{:});
            
        end
        
        [Xsh,Xwin,makeplot,iter] = align(me,Xwin,Gpart,maxiter,makeplot,compno)
        FFXpart = update_bispectrum(me,FXs,initialize)
        [Xchop,T,segment] = chop_input(me,xin,apply_window,delay,segment)
        
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
                if me(1).subspace_dim > 0
                    if size(Xwin,3) == size(me(1).projection,1)
                        filterfft = permute(squeeze(me(1).filterfft)*me(1).projection',[1 3 2]);
                    else
                        filterfft = me(1).filterfft;
                    end                    
                else
                    filterfft = me(1).filterfft;
                end
                    
                FXfilt(me(1).keepfreqs{1},:) = nansum(FXwin(me(1).keepfreqs{1},:,:).*repmat(filterfft(me(1).keepfreqs{1},:,:),1,size(X,2)),3); %#ok<*NANSUM>
                Xfilt = real(ifft(FXfilt));
                
%                 Xfilt = sum(real(ifft(FXwin.*repmat(me(1).filterfft(me(1).keepfreqs{1}),1,size(X,2)))),3);   
                if nargout > 3
                    mvXfilt = real(ifft(FXwin.*repmat(filterfft,1,size(X,2))));   
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
                FXshift = FXshift.*sgn;
            else
%                  Xin = X(:);
%                 filts = squeeze(me(1).filterfun);
                filts = me(1).filterfun;
                if me(1).subspace_dim >0
                    X = squeeze(X);
                    if size(X,2) == size(me(1).projection,1)
                        X = permute(squeeze(X)*me(1).projection,[1 3 2]); 
                    end
                end

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
%            chdim = find(size(in)>1,1,'last');
%            
%            chdim = chdim + (chdim==1); %In case the input is not in fact multivariate
           
           if me(1).subspace_dim==0 && size(in,3)~=size(me(1).filterfun,3) || me(1).subspace_dim>0 && size(in,3)~=size(me(1).projection,1)
               in = permute(in, [1 3 2]);
           end
%            if ~isvector(in) && size(in,chdim) ~= size(me(1).feature,3) && (me(1).subspace_dim== 0 || size(in,3) ~= size(me(1).projection,1))
% %                warning('MVHOS expected dimension %i for channels and %i for features, but size suggests they are reversed.\nThese will be exchanged now. In the future make sure the dimensions are correctly ordered,\nas this would have been missed if the number of features and channels happened to coincide.',chdim,chdim-1)
%                in = permute(in, [1:chdim-2 chdim chdim-1]);
%            end
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
            wf = me(1).waveform;%*me(1).projection';
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
          
            if me(1).subspace_dim > 0 && size(Xrec,2) == size(me(1).projection,2)
                Xrec = permute(squeeze(Xrec)*me(1).projection',[1 3 2]);
            end
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
        
    end

end



