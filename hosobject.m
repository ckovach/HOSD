classdef hosobject < handle
   
    % Class implementing higher-order spectral filtering based on Kovach
    % 2018.
    
    properties
       order = 3;
       freqs       
       B = 0;
       Bpart = 0;
       D = 1;
       PSD = 0;
      
       
       
       sampling_rate = 1;
       normalization = 'awplv';
       hos_learning_rate = .01;
       filter_adaptation_rate = .02;
%        learningrate = .02; % Asymptotic learning rate
       burnin = 20;
       window_number = 0;
       poverlap = .5;
      
       do_bsp_update = true;
       do_wave_update = true;
       do_filter_update = true;
       thresh = 0;
       threshtemp = 1;
       threshold_type = 'hard';
       keepfreqs
       pdonly = true;
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
    end
    
    methods
       
        function me = hosobject(order,N,sampling_rate,lowpass,freqs,varargin)
            
            warning('THIS SCRIPT IS UNDER DEVELOPMENT AND PROBABLY DOESN''T WORK RIGHT NOW')
            if nargin >1 && ~isempty(N)
                me.buffersize = N;
            end
            
            me.highpassval = 2/N;
            if nargin > 2 && ~isempty(sampling_rate)
                me.sampling_rate=sampling_rate;
                me.lowpassval = me.lowpassval*sampling_rate;
                me.glowpassval = me.glowpassval*sampling_rate;
                me.highpassval = me.highpassval*sampling_rate;
            end
            if nargin > 3 && ~isempty(lowpass)
                me.lowpassval=lowpass;
            end
             if nargin < 5 || isempty(freqs)
                freqs = fftfreq(me.bufferN)*me.sampling_rate;
            end
            if isnumeric(freqs)
                freqs = {freqs};
            end
            if nargin < 1 || isempty(order)
                order = max(length(freqs),me.order);
            end
            
            if order > length(freqs)
                freqs(end+1:order) = freqs;
            end
            me.order = order;
            me.freqs = freqs;
            me.update_frequency_indexing
            me.reset();
            me.G = ones(sum(me.keepfreqs{1}),1);
            
        end
        function reset(me)
            me.window_number = 0;
            me.waveform = zeros(me.bufferN,1);
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
            me.window_number=0;
        end
        function update_frequency_indexing(me)
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
            me.freqindx = freq2index(freqs,order,lowpass,highpass,keepfreqs,me.pdonly); %#ok<*PROP>
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
          bias = me.BIASnum./me.D.^2;
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
            
        end
        function set.filterfun(me,in)
           
            if length(in)>me.bufferN
                warning('Filer function size does not match current buffer. Filter will be truncated.')
               in(me.bufferN+1:end)=[];
            elseif length(in)<me.bufferN;
                warning('Filer function size does not match current buffer. Filter will be padded.')
                in(end+1:me.bufferN) = 0;
            end
            F =fft(ifftshift(in));
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
        
        function out = get.feature(me)
           out = ifftshift(me.waveform); 
        end
        %%%%%%%%
        function [Xfilt,FXshift] = apply_filter(me,X,varargin)
            if length(X) == me.bufferN
                FXwin = fft(repmat(me.win,1,size(X,2)).*X);
                Xfilt = real(ifft(FXwin.*repmat(me.filterfft,1,size(X,2))));   
            else
                 Xin = X;
                Xin(end+me.bufferN,:) = 0;
                Xfilt = filter(me.filterfun,1,Xin);
                Xfilt = Xfilt(floor(me.bufferN/2)+1:end-ceil(me.bufferN/2));
      
            end
            if nargout >1
                [mx,mxi] = max(Xfilt);
%                 FX = fft(X);
                delt = me.radw*me.sampt(mxi);
                FXshift = exp(1i*delt).*FXwin;
                me.delay = delt;
            end
        end
       
        %%%%%%%
        function out = get.H(me)
            
%             B = me.B(me.freqindx.remap); %#ok<*PROP>
%             B(me.freqindx.PDconj) = conj(B(me.freqindx.PDconj));
%             out = conj(B)./me.D(me.freqindx.remap).^2;
         BC = me.B./(me.D+eps);
           bias = me.BIASnum./(me.D.^2+eps);
          BC = (abs(BC)-bias).*BC./(abs(BC)+eps);
          H = BC./(me.D+eps);
          H = H(me.freqindx.remap);
          H(me.freqindx.PDconj) = conj(H(me.freqindx.PDconj));
          out = conj(H);
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
     
               %%% Remove linear phase trend so the energy of the filter
               %%% is more-or-less centered
                  dph = G(2:end).*conj(G(1:end-1));
               arg = @(x)atan2(imag(x),real(x));
               mdph = round(arg(sum(dph)./sum(abs(dph)))*length(G)/(2*pi));
               
               linphase = exp(-1i*mdph*fftfreq(length(G))'*2*pi);
              me.Bpart = me.Bpart.*[linphase(me.freqindx.Is(:,1));0];
              
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
            end
            if me.bufferPos==me.bufferN
                me.get_input(me.inputbuffer);
                me.bufferPos = 0;
            end
            getsnip = min(me.bufferN-me.bufferPos,length(snip));
            me.inputbuffer(getsnip+1:end) = me.inputbuffer(1:end-getsnip);
            me.inputbuffer(1:getsnip) = snip;
            me.bufferPos = me.bufferPos + getsnip;
        end
        
        function get_input(me,xin)
           
            % Process input if length is >= buffer size, else add to buffer.
            nxin = numel(xin);
            
            if nxin >= me.bufferN
                me.bufferPos = 0; % Discard the buffer
                stepn = round(me.poverlap*me.bufferN);
                nget = nxin - me.bufferN+1;
                tindx = (0:me.bufferN-1)';
                wint = (1:stepn:nget);
            
                T = repmat(tindx,1,length(wint))+repmat(wint,length(tindx),1);
                Xchop = xin(T);
                
                me.do_updates(Xchop)
                    
                snip = xin(T(end)+1:numel(xin));
                if ~isempty(snip)
                    me.write_buffer(snip);
                end
            else
                me.write_buffer(xin);
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
            else
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
        function Xrec = reconstruct(me,X)
            
            if nargin < 2
                Xin = me.inputbuffer;
            end
             Xfilt = me.apply_filter(X);
           
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
