classdef hosobject < handle
   
    % Class implementing higher-order spectral filtering based on Kovach
    % 2018.
    %
    % Usage:
    %        hos = hosobject(order,N,sampling_rate,lowpass,freqs,varargin)
    %
    % Inputs: 
    %       order - order (default = 3)
    %
    
    properties
       order = 3;
       freqs       
%        B = 0;
%        Bpart = 0;
%        D = 1;
       PSD = 0;
      
       
       
       sampling_rate = 1;
       normalization = 'awplv';
       hos_learning_rate = .01;
       filter_adaptation_rate = .02;
%        learningrate = .02; % Asymptotic learning rate
       burnin = 20;
       window_number = 0;
       poverlap = .5;
      
       do_update = true;
       do_bsp_update = true;
       do_wave_update = true;
       do_filter_update = true;
       thresh = 0;
       threshtemp = 1;
       threshold_type = 'hard';
       keepfreqs
       pdonly = true;
       dat = [];
      
     end
  
    properties (GetAccess = public, SetAccess=protected)
       
        inputbuffer = [];
        outputbuffer = [];
        reconbuffer = [];
        residualbuffer = [];
        shiftbuffer = [];
        thresholdbuffer=[];
        freqindx = [];
        waveform = [];
        bufferPos = 0;
        sumlr=0;
        sumlr2=0;
        radw = [];
        sampt = [];
        delay = 0;
    end
    
    properties (Access = private)
      bufferN = 1024;
      G = [];

      wintype = 'hann';
      win = hann(1024);
      BCpart = 0;
       highpassval = 0;
       lowpassval = .5; %%% Lowpass on the edges (max freq)
       glowpassval = .5; %%% Global lowpass
      BIASnum = 0;
        
       Bval = 0;
       Bpartval = 0;
       Dval = 1;
    end
    properties (Dependent = true)
        
%         inputfft
%         inputfiltered  
%         input_shifted
        buffersize
        filterfun
        filterfft
%         current_learning_rate;
        window 
        bicoh
        H
        EDF
        highpass  ;
        lowpass ; %%% Lowpass on the edges (max freq)
        glowpass ; %%% Global lowpass
        Bfull
        BIAS
        fullmap
        feature
        current_threshold % Current adaptive threshold level
        B ;
        Bpart ;
        D ;
    end
    
    methods
       
        function me = hosobject(order,varargin)
            
            warning('THIS SCRIPT IS UNDER DEVELOPMENT AND PROBABLY DOESN''T WORK RIGHT NOW')
            if nargin < 1 
                return
            elseif nargin == 1
                me.order = order;
                return
            else
                me.order = order;
            end
                
            me.initialize(order,varargin{:});
        end
        
      function initialize(me,N,sampling_rate,lowpass,freqs,freqindex,varargin)
      
            if ~isscalar(N)
                X = N;
                N = size(X,1);
            else
                X = [];
             end
             if nargin >1 && ~isempty(N)
                me(1).buffersize = N;
             end
            if nargin < 6
                freqindex = [];
            end
            me(1).highpassval = 2/N;
            if nargin > 2 && ~isempty(sampling_rate)
                me(1).sampling_rate=sampling_rate;
                me(1).lowpassval = me(1).lowpassval*sampling_rate;
                me(1).glowpassval = me(1).glowpassval*sampling_rate;
                me(1).highpassval = me(1).highpassval*sampling_rate;
            end
            if nargin > 3 && ~isempty(lowpass)
                me(1).lowpassval=lowpass;
            end
             if nargin < 5 || isempty(freqs)
                freqs = fftfreq(me(1).bufferN)*me(1).sampling_rate;
            end
            if isnumeric(freqs)
                freqs = {freqs};
            end
            if length(freqs) > me(1).order
                me(1).order = length(freqs);
            end
            
            if me(1).order > length(freqs)
                freqs(end+1:me(1).order) = freqs;
            end
%             me.order = order;
            me(1).freqs = freqs;
            me(1).update_frequency_indexing(freqindex)
            me(1).reset();
            me(1).G = ones(sum(me(1).keepfreqs{1}),1);
            
            if length(me)>1
                me(2:end).initialize(N,sampling_rate,lowpass,freqs,freqindex,varargin{:});
            end
            
            if ~isempty(X)
                me.get_block(X);
            end
        end
        function reset(me)
            me.window_number = 0;
            me.sumlr =0;
            me.sumlr2 = 0;
            z=zeros(me.bufferN,1);
            me.inputbuffer = z;
            me.outputbuffer = z;
            me.waveform = z;
            me.shiftbuffer = z;
            me.thresholdbuffer=z;
%             me.G = ones(size(z));
            me.bufferPos = 0;
            me.B(:)=0;
            me.G(:)=1;
            me.window_number=0;
        end
        function update_frequency_indexing(me,freqindex)
            lowpass = me.lowpassval;
            order = me.order;
            freqs=me.freqs;
            if length(lowpass)< order-1
                lowpass(end+1:order-1) = lowpass(end);
            end
            if length(lowpass)< order
                lowpass(order) = me.glowpassval;
            end
            
            highpass = me.highpassval;
            if length(highpass)< order
                highpass(end+1:order) = highpass;
            end
            
            keepfreqs={};
            for k = 1:length(me.freqs)                
                keepfreqs{k} =(abs(me.freqs{k})<=lowpass(k)&abs(me.freqs{k})>highpass(k));                 %#ok<*AGROW>
%                 freqs{k} = freqs{k}(keepfreqs{k});
               % freqindex{k} = find(keepfreqs{k});
            end    
            me.keepfreqs = keepfreqs;
            %%% Initialize the indexing   
            if nargin < 2 || isempty(freqindex)
                freqindx = freq2index(freqs,order,lowpass,highpass,keepfreqs,me.pdonly); %#ok<*PROPLC,*PROP>
            end
            
            me.freqindx  = freqindx;
                
            Z =zeros(size(me.freqindx.Is,1)+1,1);
            me.B = Z; 
            me.Bpart = Z;
            me.D = Z;
             me.BIASnum=Z;
            me.G= ones(sum(me.keepfreqs{1}),1);
             me.reset;
        end
        function [lradj,lr] = learningfunction(me,learningrate,m,burnin)
            if nargin < 3 || isempty(m)
                m = 1;
            end
            if nargin < 4 || isempty(burnin)
                burnin = me.burnin;
            end
            lr = exp(-me.window_number./burnin)./(me.window_number+1) + (1-exp(-me.window_number./burnin))*learningrate;
            lradj = (1-(1-lr)^m);

        end
        
        function out = get.Bfull(me)
          out = me.B(me.freqindx.remap);
          out(me.freqindx.PDconj) = conj(out(me.freqindx.PDconj));
        end
        function BC = get.bicoh(me)
            
          BC = me.B./me.D;
           bias = me.BIASnum./(me.D.^2+eps);
          BC = (abs(BC)-bias).*BC./(abs(BC)+eps);
          BC = BC(me.freqindx.remap);
          BC(me.freqindx.PDconj) = conj(BC(me.freqindx.PDconj));
          
        end
        function out = get.BIAS(me)        
          bias = me.BIASnum./(me.D.^2+eps);
          out = bias(me.freqindx.remap);
        end
        function out = get.fullmap(me)
           out = me.freqindx.remap; 
        end
        function set.lowpass(me,a)
           me.lowpassval = a;
            me.update_frequency_indexing;
        end
        function set.glowpass(me,a)
           me.glowpassval = a;
            me.update_frequency_indexing;
        end
        function set.highpass(me,a)
           me.highpassval = a;
            me.update_frequency_indexing;
        end
        function out = get.lowpass(me)
          out = me.lowpassval ;
        end
        function out = get.glowpass(me)
           out = me.glowpassval;
        end
        function out = get.highpass(me)
          out =  me.highpassval;
        end
        function out = get.buffersize(me)
            out = me.bufferN;
        end
        function set.buffersize(me,N)
           me.bufferN = N;
           me.win = window(me.wintype,N); 
           me.radw = ifftshift((0:me.bufferN-1)' - floor(me.bufferN/2))/me.bufferN*2*pi;
           me.sampt = ifftshift((0:me.bufferN-1) - floor(me.bufferN/2)'); 
           me.reset;
        end
        function out = get.window(me)
           out =  me.wintype;
        end
        function set.window(me,win)
           me.wintype=win;
            me.win = window(win,me.bufferN); %#ok<*CPROPLC>
        end
        function out = get.filterfft(me)
           
            out = zeros(me.bufferN,1);
            out(me.keepfreqs{1}) = me.G;
            
        end
        function out = get.filterfun(me)
           
            F = me.filterfft;
            out = ifftshift(real(ifft(F)));
%             out = real(ifft(F));
            
        end
        function set.filterfun(me,in)
           
            if length(in)>me.bufferN
                warning('Filer function size does not match current buffer. Filter will be truncated.')
               in(me.bufferN+1:end)=[];
            elseif length(in)<me.bufferN;
                warning('Filer function size does not match current buffer. Filter will be padded.')
                in(end+1:me.bufferN) = 0;
            end
%             F =fft(fftshift(in));
            F =fft((in));
            me.filterfft = F;
            
        end
        function set.filterfft(me,in)
           
            me.G = in(me.keepfreqs{1}) ;
            
        end
%         %%%%%%%
%         function out = get.current_learning_rate(me)
%             out = me.learningfunction;
%         end
        %%%%%%%%
        function out = get.EDF(me)
            
            % Cumulative effective degrees of freedom based on learning rate
            out = me.sumlr./me.sumlr2;
        end
        %%%%%%%
        function set.EDF(me,in)
            
            % When setting the EDF, take it to be a simple sample size
            me.sumlr = 1;
            me.sumlr2 = 1./in;
            me.window_number = in;
        end
        %%%%%%%
        function out = get.feature(me)
           out = ifftshift(me.waveform); 
        end
         %%%%%%%
        function set.feature(me,in)
           me.waveform = fftshift(in); 
        end
        %%%%%%%%
        function [Xfilt,FXshift] = apply_filter(me,X,apply_window,varargin)
            if nargin<3
                apply_window = true;
            end
            if length(X) == me.bufferN
                if apply_window
                    win = me.win;
                else
                    win =1;
                end
                Xwin = fftshift(repmat(win,1,size(X,2)).*X,1);
                FXwin = fft(Xwin);
%                 FXwin = fft(X)';
                Xfilt = real(ifft(FXwin.*repmat(me.filterfft,1,size(X,2))));   
            else
                 Xin = X;
                Xin(end+me.bufferN,:) = 0;
                Xfilt = filter(me.filterfun,1,Xin);
                Xfilt = Xfilt(ceil(me.bufferN/2)+1:end-floor(me.bufferN/2));
      
            end
            if nargout >1
                [mx,mxi] = max(Xfilt);
%                 FX = fft(X);
                samptc=(me.sampt);
                 dt = samptc(mxi);
                delt = me.radw*dt;
                FXshift = exp(1i*delt).*FXwin;
                me.delay = dt;
            end
        end
       
        function out = get.B(me)
           out = me.Bval; 
        end
        function out = get.D(me)
           out = me.Dval; 
        end
        function out = get.Bpart(me)
           out = me.Bpartval; 
        end
        function set.B(me,in)
            if min(size(in))==1
                me.Bval = in; 
            else
                me.Bval = [in(me.freqindx.reduce);0];
            end
        end
        function set.D(me,in)
            if min(size(in))==1
              	me.Dval=in; 
            else
                me.Dval  = [in(me.freqindx.reduce);0];
            end
        end
        function  set.Bpart(me,in)
            if min(size(in))==1
               	me.Bpartval=in; 
            else
               me.Bpartval = [in(me.freqindx.reduce);0];
            end
        end
        %%%%%%%
        function out = get.H(me)
            
%             B = me.B(me.freqindx.remap); %#ok<*PROP>
%             B(me.freqindx.PDconj) = conj(B(me.freqindx.PDconj));
%             out = conj(B)./me.D(me.freqindx.remap).^2;
         BC = me.B./(me.D+eps);
           bias = me.BIASnum./(me.D.^2+eps);
           bias(isnan(bias))=0;
          BC = (abs(BC)-bias).*BC./(abs(BC)+eps);
          H = BC./(me.D+eps);
          H = H(me.freqindx.remap);
          H(me.freqindx.PDconj) = conj(H(me.freqindx.PDconj));
          out = conj(H);
        end
        %%%%%%%
        function out = xfilt(me,in)
           if nargin < 2
               in = me.dat;
           end
           out = me(1).apply_filter(in);  
           if length(me)>1
               out = [out,me(2:end).xfilt(in-me(1).xrec(in))];
           end
        end
        %%%%%%%
        function out = ximp(me,in)
           % Just get the thresholded times (for backward compatibility)
           if nargin < 2
               in = me.dat;
           end
           out = me.xthresh(in)>0;  
%            if length(me)>1
%                out = [out,me(2:end).ximp(in-me(1).xrec(in))];
%            end         

        end
        function out = xthresh(me,in)
             % Get the thresholded data 
           if nargin < 2
               in = me.dat;
           end
           out = sparse(me(1).filter_threshold(me(1).apply_filter(in)));  
           if length(me)>1
               out= [out,me(2:end).xthresh(in-me(1).xrec(in))];
           end
         

        end
        %%%
        function out = xrec(me,in)
           if nargin < 2
               in = me.dat;
           end
           
           out = me(1).reconstruct(in); 
           if length(me)>1
               out = [out,me(2:end).xrec(in-out(:,1))];
           end
        end
        %%%%%%%
        function update_bispectrum(me,FX)
            
            %Write now this updates in chunks with temporal decay weighting
            %applied only serially. That is, a simple average is obtained
            %for each chunk, which is then 
            m = size(FX,2);

            if isempty(FX)
                return
            end
%             Xwin = Xin.*me.win;
%             FX = fft(Xwin);
           % FFX = 1;
            FFX = conj(FX(me.freqindx.Is(:,me.order),:));
            for k = me.order-1:-1:1
                if k == 1
                    FFXpart = FFX;
                end
                FFX = FFX.*FX(me.freqindx.Is(:,k),:);
            end
            
            BX = mean(FFX,2);
            XPSD = mean(abs(FX).^2,2);
            BXpart = mean(FFXpart,2);
            BX(end+1,1) = 0;
            XPSD(end+1,:) =0;
            BXpart(end+1,:) = 0;
            
            [lradj,lr] = me.learningfunction(me.hos_learning_rate,m);
             fflr = me.learningfunction(me.filter_adaptation_rate,m,1./me.filter_adaptation_rate);
%             fflr = (1-(1-me.filter_adaptation_rate)^m);
            
            %%% Adjust the learning rate according to the number of samples
            %%% in Xin, giving the sample weight the same total as if it
            %%% had been added serially.
             asympedf = 2./lr-1; %Asymptotic EDF
       
            lrbias = 1./asympedf*(1-(1-lr).^(2*m)); % "Learning rate" for the sum of squared weights in the bias term
        
            me.B = (1-lradj)*me.B + lradj*BX;
            me.Bpart = (1-fflr)*me.Bpart + fflr*BXpart;
            me.PSD = (1-lradj)*me.PSD + lradj*XPSD;
            me.sumlr = me.sumlr*(1-lradj) + lradj;
            me.sumlr2 = me.sumlr2*(1-lr).^(2*m) + lrbias;
%            me.sumlr2 = (me.sumlr2-1./asympedf)*(1-me.current_learning_rate).^(2*m) + lrbias;
            
            switch me.normalization
                case 'awplv'
                    NX = mean(abs(FFX),2);
                    NX(end+1,1) = 0;
                    XbiasNum = sum(abs(FFX).^2,2)./m^2;
                    XbiasNum(end+1,1) = 0;
                    me.BIASnum = me.BIASnum.*(1-lr).^(2*m) + lrbias*XbiasNum;
                    me.D = (1-lradj)*me.D + lradj*NX+eps;
                   
                case {'bicoh','bicoherence'}
                    XBCpart = mean(FFXpart,2);
                    XBCpart(end+1,1) = 0;
                    me.BCpart = (1-lradj)*me.BCpart + lradj*XBCpart;                    
                    me.D = sqrt(me.BCpart.*me.PSD(me.freqindx.Is(:,1)))+eps;
            end
            
         
        end
        
        function update_filter(me)
               
               Bpart = me.Bpart(me.freqindx.remap);
               Bpart(me.freqindx.PDconj) = conj(Bpart(me.freqindx.PDconj));
               G = sum(Bpart(:,:).*me.H(:,:),2);
     
%                %%% Remove linear phase trend so the energy of the filter
%                %%% is more-or-less centered
%                   dph = G(2:end).*conj(G(1:end-1))./(abs(G(1:end-1))+abs(G(2:end))+eps)*2;
%                arg = @(x)atan2(imag(x),real(x));
%                mdph = (arg(sum(dph)./sum(abs(dph)))*length(G)/(2*pi));
%                
%                linphase = exp(-1i*mdph*fftfreq(length(G))'*2*pi);
%               me.Bpart = me.Bpart.*[linphase(me.freqindx.Is(:,1));0];
              
          %    me.G = Gshift(me.keepfreqs{1}(abs(me.freqs{1})<=me.lowpass(1)));
                me.G = G(me.keepfreqs{1}(abs(me.freqs{1})<=me.lowpass(1)));
        end
        
        function update_waveform(me,Xsh)
           m = size(Xsh,2);
           lradj = me.learningfunction(me.filter_adaptation_rate,m,1./me.filter_adaptation_rate);
%             lradj = (1-(1-me.filter_adaptation_rate)^m);
           me.waveform = me.waveform*(1-lradj) + mean(Xsh,2)*lradj; 
        end
        
        function write_buffer(me,snip)
            if length(snip)>me.bufferN
                me.get_input(snip);
                return
            end
            if me.bufferPos==me.bufferN
                me.get_input(me.inputbuffer);
                me.bufferPos = 0;
                return
            end
            getsnip = min(me.bufferN-me.bufferPos,length(snip));
            me.inputbuffer(length(snip)+(1:end-length(snip))) = me.inputbuffer(1:me.bufferN-length(snip));
            me.inputbuffer(1:length(snip)) = snip;
            me.bufferPos = me.bufferPos + length(snip);
        end
        
        function get_input(me,xin)
           
            % Process input if length is >= buffer size, else add to buffer.
            nxin = numel(xin);
            
            if nxin >= me(1).bufferN
                me(1).bufferPos = 0; % Discard the buffer
                if  size(xin,1) ~=me(1).bufferN 
                    stepn = round(me(1).poverlap*me(1).bufferN);
                    nget = nxin - me(1).bufferN+1;
                    tindx = (0:me(1).bufferN-1)';
                    wint = (1:stepn:nget);

                    T = repmat(tindx,1,length(wint))+repmat(wint,length(tindx),1);
                    Xchop = xin(T);
                                
                    snip = xin(T(end)+1:numel(xin));
                    if ~isempty(snip)
                        me(1).write_buffer(snip);
                    end
                else
                    Xchop = xin;
                end
                me(1).do_updates(Xchop)
        
            else
                me(1).write_buffer(xin);
            end    
           
            if length(me)>1
               xrec = me(1).reconstruct(xin);
               me(2:end).get_input(xin-xrec);
            end
            
            
            
        end
        
        function get_block(me,xin,maxiter,makeplot,compno)
           
            % Fit a block of data all at once
            % Process input if length is >= buffer size, else add to buffer.
            if nargin < 5
                compno = 1;
            end
            if nargin < 4
                makeplot = true;
            end
            if nargin < 3 || isempty(maxiter)
                maxiter = 100;
            end
            nxin = numel(xin);
            
            if nxin >= me(1).bufferN
                me(1).bufferPos = 0; % Discard the buffer
                if  size(xin,1) ~=me(1).bufferN 
                    stepn = round(me(1).poverlap*me(1).bufferN);
                    nget = nxin - me(1).bufferN+1;
                    tindx = (0:me(1).bufferN-1)';
                    wint = (1:stepn:nget);

                    T = repmat(tindx,1,length(wint))+repmat(wint,length(tindx),1);
                    Xchop = xin(T);
                                
                    snip = xin(T(end)+1:numel(xin));
                    if ~isempty(snip)
                        me(1).write_buffer(snip);
                    end
                else
                    Xchop = xin;
                end
                del = Inf;
                tol =1; % Stop when the average shift is less than 1 sample
                k = 0;
                olddt2 = 0;
                Xsh = Xchop.*repmat(me(1).win,1,size(Xchop,2));
                fprintf('\nComponent %3i Iter %3i',compno,0)
                while del >tol && k < maxiter                    
                    if ishandle(makeplot)
                        set(makeplot,'cdata',Xsh);
                        title(sprintf('Component %3i, Iter. %3i, Mean shift = %0.2fs',compno,k,del/me(1).sampling_rate));
                        drawnow
                    elseif islogical(makeplot) && makeplot
                        figure,
                        makeplot = imagesc([],fftshift(me(1).sampt)/me(1).sampling_rate,Xsh);
                        title(sprintf('Component %03i, Iter. %03i, Mean shift = %0.2fs',compno,k,del/me(1).sampling_rate));
                        drawnow
                    end
                    
                    k=k+1;
                    fprintf('\b\b\b%03i',compno,k)
                    
                    olddt = me(1).delay;
                    me(1).get_input(Xsh)
                    [~,FXsh] = me(1).apply_filter(Xsh,false);
                    Xsh = ifftshift(ifft(FXsh),1);
                    newdt = me(1).delay;
                    % checks two and one step back to reduce getting
                    % trapped at points of cyclical stability.
                    del = min(sqrt(mean((olddt-newdt).^2)),sqrt(mean((olddt2-newdt).^2)));
                    olddt2 = olddt;
                end
            else
                me(1).write_buffer(xin);
            end    
           
            if length(me)>1
               xrec = me(1).reconstruct(xin);
               me(2:end).get_block(xin-xrec,maxiter,makeplot,compno+1);
            end
            
            
            
        end
        
        
        
        function [Xthresh,Xcs,trialthresh] = filter_threshold(me,Xfilt,thresh)
            
            % Apply a moment-based threshold
            if nargin < 3
                thresh= me.thresh;
            end
            
             Xcent = zscore(Xfilt);
            %Xcent = Xfilt;
            Xmom = Xcent.^me.order;
            
            if mod(me.order,2)==0
                Xmom = Xmom - me.order*repmat(mean(Xcent.^2).^(me.order/2),size(Xmom,1),1);
            end
            if size(Xfilt,1) == me.bufferN;
                 trialthresh = me.current_threshold;
            elseif me.order ==3
                % For the bispectrum compute normalized skewness
%                 keepsamples = ones(size(Xcent));
                  srt = sort(Xcent);
                outlier_threshold = 5;
               keepsamples = ~isnan(iterz(srt,outlier_threshold,-1)); % Suppress extreme negative outliers           
                m1 = cumsum(srt.*keepsamples)./cumsum(keepsamples); % cumulative mean on sorted peaks
                m2 = cumsum(srt.^2.*keepsamples)./cumsum(keepsamples); % cumulative 2nd moment
                m3 = cumsum(srt.^3.*keepsamples)./cumsum(keepsamples); % cumulative 3rd moment
                %  Third cumulant
                c3 = m3 - 3*m2.*m1 + 2*m1.^3; % Third cumulant on sorted peaks
          
                keepsrt = srt>0 & c3> me.thresh ;
                detect = any(keepsrt);
                trialthresh = sum ((diff(keepsrt)>0).*srt(2:end,:)).^me.order;
                trialthresh(~detect) = Inf;
            else
                % For now, apply a simple threshold on the moment for
                % orders > 3. This should be improved to use the proper
                % cumulant.
                Xsrt = sort(Xmom);
                Xcs = cumsum(Xsrt)>thresh;
                detect =any(Xcs);
                trialthresh = sum(Xsrt(2:end,:).*diff(Xcs));
                trialthresh(~detect) = Inf;
            end
            switch me.threshold_type
                case 'hard'
                    THR = (Xmom>=repmat(trialthresh,size(Xmom,1),1));
                case 'soft'
                    sigm = @(x)1./(1+exp(-x));
                    THR = sigm(me.threshtemp*(Xmom-repmat(trialthresh,size(Xmom,1),1)));
                otherwise
                    error('Unrecognized threshold type')
            end
            Xthresh = Xfilt.*THR;
            
        end
        function [Xrec,Xfilt] = reconstruct(me,X)
            
            if nargin < 2
                Xin = me.inputbuffer;
            end
            Xfilt = me.apply_filter(X);
            if length(X) == me.bufferN
                Xfilt = ifftshift(Xfilt,1);
            end
            
           % Xfilt = Xfilt(floor(me.bufferN/2):end-ceil(me.bufferN/2));
            FXthresh =fft(me.filter_threshold(Xfilt));
            wf = ifftshift(me.waveform);
            wf(size(FXthresh,1)) = 0;
            wf = circshift(wf,-floor(me.bufferN/2));
            featfft = fft(wf);
            Xrec = real(ifft(FXthresh.*repmat(featfft,1,size(X,2))));
%             Xrec = Xrec(floor(me.bufferN/2)+1:end-ceil(me.bufferN/2));
           a= sum(abs(Xrec(:)).^2);
           if a > 0
             Xrec = Xrec*(X(:)'*Xrec(:))./a; % Scale to minimize total mse.
           end
           if nargin < 2
                me.reconbuffer = Xrec;
            end
        end
        function do_updates(me,X)
            
            if ~me.do_update
               warning('Updating is currently disabled. Set do_update = true to enable.') 
               return
            end
            [Xfilt,FXsh] = me.apply_filter(X);
            getwin = me.update_criteria(Xfilt);
            if isempty(getwin)
                return
            end
            if me.do_bsp_update
               me.update_bispectrum(FXsh(:,getwin)); 
            end
            if me.do_filter_update
               me.update_filter(); 
               Xsrt = mean(sort(zscore(Xfilt(:,getwin)).^me.order),2);
               lradj = me.learningfunction(me.filter_adaptation_rate,sum(getwin));
               me.thresholdbuffer = me.thresholdbuffer*(1-lradj) + lradj*cumsum(Xsrt);
            end
            if me.do_wave_update
                Xsh = real(ifft(FXsh(:,getwin)));
                me.update_waveform(Xsh); 
            end
            me.window_number = me.window_number+sum(getwin);
            me.outputbuffer = mean(Xfilt,2);
            me.shiftbuffer = real(ifft(mean(FXsh,2)));
           
        end
        
        function out = update_criteria(me,Xfilt)
            out = true(1,size(Xfilt,2)); % Placeholder for now
        end
        
        function out =get.current_threshold(me)
           xsrt = diff(me.thresholdbuffer);
           threshi = find(diff(me.thresholdbuffer>me.thresh));
           if isempty(threshi)
               threshi=length(xsrt);
           end
           out = xsrt(threshi);
        end
    end
    
end

function out = fftfreq(N)
    out = ifftshift((0:N-1)-floor(N/2))/N;
end
