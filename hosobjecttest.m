
classdef hosobjecttest < handle
   
    % Class implementing higher-order spectral filtering based on Kovach
    % and Howard 2019.
    %
    % Usage: 
    %   To create a hos object
    %        hos = hosobject(order,N,sampling_rate,lowpass)
    %
    %   To fit the object to a block of data (offline mode)
    %        hos.get_block(data, [maxiter=50])
    %
    %   To add a segment of data in computing a running average (online mode):
    %        hos.get_input(data)
    %   
    %   To initialize an M component decomposition:
    %        hos(M) = hosobject;
    %        hos.initialize(N,sampling_rate,lowpass)
    %   
    % Inputs: 
    %       HOS order - order (default = 3)
    %       N - buffer length in samples used to compute HOS
    %       sampling_rate - sample rate
    %       lowpass - lowpass cutoff
    %       data  -  input data in the form of samples x segments. If data
    %               is a single column vector it will be segmented into
    %               overlapping N point segments.
    %       maxiter - maximum iterations (default - 25)
    %
    % Outputs: 
    %       hos.waveform - Recovered feature waveform(s)
    %       hos.filterfun - Feature detection filter(s)
    %       hos.bicoh  -  Bicoherence of the input signal (or polycoherence for orders > 3)
    %       xfilt = hos.apply_filter(data) - apply the detection filter to the data
    %       xthresh = hos.xthresh(data) - Thresholded signal used in the reconstruction.
    %       ximp = hos.ximp(data) - Suprathreshold samples (samples at which the feature is detected).
    %       xrec = hos.xrec(data) - Reconstructs the signal(s) associated with one or more features.       
    %                  
    %
    %
    % Copyright Christopher K. Kovach, University of Iowa 2018
    
    properties
       order = 3;
       freqs       
%        B = 0;
%        Bpart = 0;
%        D = 1;
       PSD = 0;
      
       
       
       sampling_rate = 1; % Sampling rate of the input data
       normalization = 'awplv'; % Normalization used in computing bi-(or poly-)coherence
       hos_learning_rate = .01; % Learning rate for online mode
       filter_adaptation_rate = .02; % Filter adaptation rate for online mode
%        learningrate = .02; % Asymptotic learning rate
       burnin = 20;
       window_number = 0;  % Nunmber of window processed
       poverlap = .5;      % Default overlap between adjacent windows as the 
                           % proportion of step size to window width (poverlap = 1 
                           % is no overlap, poverlap = 2 interleaves a gap of 1
                           % window duration between windows.)
      
       do_update = true;
       do_bsp_update = true;
       do_wave_update = true;
       do_filter_update = true;
       adjust_lag = true; % Automatically apply a circular shift to the filter and waveforms to center the energy in both
       lag = 1; % A phasor representing the amount of circularshift added to the filter estimate (1 = no shift, +/-1i = max shift)
       thresh = 0;
       threshtemp = 1;
       threshold_type = 'hard';
       threshold_order = [];
       keepfreqs
       pdonly = true;
       dat = [];
       avg_delay = 1; % Average delay is stored as a phasor because averaging is in the circular domain.
       outlier_threshold = 5;
       fftN = 1024;
       regstat = [];
       regweight=[];
       sampweight = [];
       Imats = {};
       Iconjmats = {};
       segment = struct('wint',[],'wintadj',[],'Trange',[],'fs',1,'discarded',[]);
       
     end
  
    properties (GetAccess = public, SetAccess=protected)
       
        inputbuffer = []; % Current input buffer
        outputbuffer = []; % Current outputbuffer
        reconbuffer = [];
        residualbuffer = [];
        shiftbuffer = [];
        thresholdbuffer=[];
        freqindx = [];
        bufferPos = 0;
        sumlr=0;
        sumlr2=0;
        radw = [];  % Vector of sample frequencies in radian units
        sampt = []; % Vector of sample indices
        delay = 0;
        waveftlag = [];
      
        win = sasaki(1024);
      G = []; 

    end
    
    properties (Access = protected)
      bufferN = 1024;
      

%      wintype = 'hann'; % Default window type
      wintype = 'sasaki'; % Default window type
 
      BCpart = 0;
       highpassval = 0;
       lowpassval = .5; %%% Lowpass on the edges (max freq)
       glowpassval = .5; %%% Global lowpass
       slowpassval = .5; %%% freq. pair lowpass
       shighpassval = 0; %%% freq. pair highpass
       xlowpassval = .5; %%% OR lowpass
       xhighpassval = 0; %%% OR lowpass
      BIASnum = 0;
        
       Bval = 0;
       Bpartval = {};
       Dval = 1;
       padN = 0;       
       regval=[];
       
       use_partial_delay_method = true;
            %%% If false, uses the older method which computed HOS  before
            %%% filer estimation, which is unnecessary in iterated estimation.
            %%% If true, then the partial_delay_filters are estimated for each
            %%% segment and the detection filter is obtained by shifting and 
            %%% averaging the segment filters. This option is included for the 
            %%% sake of debugging and will be removed eventually. 
        
       do_indexing_update = true;
    end
    properties (Dependent = true)
        
%         inputfft
%         inputfiltered  
%         input_shifted
        buffersize
        waveform; %% The feature waveform
        wavefft;
        filterfun %% The ifftshifted filter function 
        filterfft %% FFT of the filterfunction
%         current_learning_rate;
        filterftlag;
        window %% Window used prior to calculating estimates
        bicoh
        partialbicoh
        H
        EDF
        highpass  ;
        lowpass ; %%% Lowpass on the edges (max freq)
        glowpass ; %%% Global lowpass
        shighpass ; %%% Highpass applied to sums of frequency pairs for order > 3.
        slowpass ; %%% Lowpass applied to sums of frequency pairs for order > 3.
        
        xlowpass ; %%% "OR" lowpass and highpass: include regions in which ANY of the frequencies meet the criterion 
        xhighpass ;%%% This is useful to design filters selective for the interior or exterior regions of
                   %%% the bispectrum.
        Bfull
        BIAS
        fullmap
        feature
        current_threshold % Current adaptive threshold level
        B ;
        Bpart ;
        D ;
        pad; %
        regressor;
    end
    
    methods
       
        function me = hosobjecttest(order,varargin)
            
            if nargin ==0
                return
            end
            if isa(order,mfilename) || isa(order,'hosminimal')|| isa(order,'hosobject')
               obj = order;
               fns = [{'BIASnum'};setdiff(properties(obj),{'BIAS','Bfull','H','bicoh','current_threshold','sampling_rate','freqindx','buffersize','filterftlag','fullmap','partialbicoh','filterfft','filterfun'})];
               
               
               if length(obj)==1
                   obj(2:length(me)) = obj;
               elseif length(obj)>length(me)
                   me(length(obj)) = hosobject;
               end
               
               
               me(1).order = obj(1).order;               
               me.initialize(obj(1).bufferN,obj(1).sampling_rate,obj(1).lowpass,obj(1).freqs,obj(1).freqindx,varargin{:})
               
               me(1).do_indexing_update = false;
               for k = 1:length(fns)  
                   if isprop(obj(1),fns{k})
                       me(1).(fns{k}) = obj(1).(fns{k});
                   end
               end
               me(1).do_indexing_update = true;
               if length(me)>1
                   me(2:end) = hosobject(obj(2:end));
               end
               return
            end
            if nargin < 1 
                return
            elseif nargin == 1
                me.order = order;
                return
            else
                me.order = order;
            end
                
            me.initialize(varargin{:});
        end
        
      function initialize(me,N,sampling_rate,lowpass,freqs,freqindex,varargin)
      
          
            if ~isscalar(N)
                X = N;
                N = size(X,1);
                me(1).do_update = true;
            else
                X = [];
             end
             if nargin >1 && ~isempty(N)
                me(1).bufferN = N;
                me(1).fftN = N;
             end
            if nargin < 6
                freqindex = [];
            end
            me(1).highpassval = 2/N;
            me(1).shighpassval = 2/N;
            if nargin > 2 && ~isempty(sampling_rate)
                me(1).sampling_rate=sampling_rate;
%                 me(1).lowpassval = me(1).lowpassval*sampling_rate;
%                 me(1).glowpassval = me(1).glowpassval*sampling_rate;
%                 me(1).highpassval = me(1).highpassval*sampling_rate;
            else
                sampling_rate = me(1).sampling_rate;
            end
            if nargin > 3 && ~isempty(lowpass)
                me(1).lowpassval=lowpass./me(1).sampling_rate;
            else
                lowpass = me(1).lowpass;
            end
             if nargin < 5 || isempty(freqs)
                freqs = fftfreq(me(1).fftN)*me(1).sampling_rate;
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
%            me(1).G = ones(sum(me(1).keepfreqs{1}),1);
            me(1).do_indexing_update = false;
            k = 1;
            me(1).do_indexing_update = false;
            while k < length(varargin)    
                me(1).(varargin{k}) = varargin{k+1};
                k=k+2;
            end
            me(1).do_indexing_update = true;
            me(1).update_frequency_indexing(freqindex)
            me(1).reset();
            
            if length(me)>1
                me(2:end).initialize(N,sampling_rate,lowpass,freqs,freqindex,varargin{:});
            end
            
            if ~isempty(X)
                me.get_block(X);
            end
        end
        function reset(me)
            
            me(1).radw = ifftshift((0:me(1).fftN - 1 )' - floor((me(1).fftN)/2))/(me(1).fftN)*2*pi;
            me(1).sampt = ifftshift((0:me(1).fftN - 1 ) - floor((me(1).fftN)/2)'); 
            me(1).window_number = 0;
            me(1).sumlr =0;
            me(1).sumlr2 = 0;
            z=zeros(me(1).bufferN,1);
            z2=zeros(me(1).fftN,1);
            me(1).inputbuffer = z;
            me(1).outputbuffer = z;
            me(1).shiftbuffer = z2;
            me(1).thresholdbuffer=z;
            me(1).PSD = [z;0];
%             me(1).G = ones(size(z));
            me(1).bufferPos = 0;
            me(1).B(:)=0;
            me(1).G(:)=1;
            me(1).D(:)=1;          
            me(1).window_number=0;
            me(1).avg_delay = 1;
            me(1).lag=1;
            me(1).Imats = {};
            me(1).Iconjmats = {};
            me(1).waveftlag=z;
            
            me(1).win = window(me(1).window,me(1).fftN);
            if me(1).fftN< me(1).bufferN || (~isempty(me(1).keepfreqs) && length(me(1).keepfreqs{1})~=me(1).bufferN)
                me(1).buffersize = me(1).bufferN;
%                 me(1).fftN = me(1).bufferN;
            end
            for k = 1:length(me(1).Bpart)
                me(1).Bpart{k}(:) = 0;
            end
      
            me(1).waveform = z2;

            if length(me)>1
                me(2:end).reset();
            end
        end
        function update_frequency_indexing(me,freqindx,mask)
           
            if ~me(1).do_indexing_update
                return
            end
            order = me.order;
            freqs=me.freqs;
            
            lowpass = me.lowpassval*me.sampling_rate;
            if length(lowpass)< order-1
                lowpass(end+1:order-1) = lowpass(end);
            end
            if length(lowpass)< order
                lowpass(order) = me.glowpassval*me.sampling_rate;
            end
            xlowpass = me.xlowpassval*me.sampling_rate;
            if length(xlowpass)< order
                xlowpass(end+1:order) = xlowpass(end);
            end
            slowpass = me.slowpassval*me.sampling_rate;
            shighpass = me.shighpassval*me.sampling_rate;
            
            xhighpass = me.xhighpassval*me.sampling_rate;
            if length(xhighpass)< order
                xhighpass(end+1:order) = xhighpass(end);
            end
            if nargin < 3 || isempty(mask)
                mask = true;
            end
            highpass = me.highpassval*me.sampling_rate;
            if length(highpass)< order
                highpass(end+1:order) = highpass(1);
            end
            
            keepfreqs={};
            for k = 1:length(me.freqs)                
                keepfreqs{k} =(abs(me.freqs{k})<=lowpass(k)&abs(me.freqs{k})>highpass(k));                 %#ok<*AGROW>
%                 freqs{k} = freqs{k}(keepfreqs{k});
               % freqindex{k} = find(keepfreqs{k});
            end    
            me.keepfreqs = keepfreqs;
            %%% Initialize the indexing   
            if nargin < 2 || isempty(freqindx)
                freqindx = freq2index(freqs,order,lowpass,highpass,keepfreqs,me.pdonly,[],mask,xlowpass,xhighpass,slowpass,shighpass); %#ok<*PROPLC,*PROP>
            end
            
            me.freqindx  = freqindx;
                
            Z =zeros(size(me.freqindx.Is,1)+1,1);
            me.B = Z; 
            me.Bpart = {};
            me.Bpart(1:me.order) = {Z};
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
           bias = sqrt(me.BIASnum./(me.D.^2+eps));
          BC = (abs(BC)-bias).*BC./(abs(BC)+eps);
          BC = BC(me.freqindx.remap);
          BC(me.freqindx.PDconj) = conj(BC(me.freqindx.PDconj));
          
        end
         function pBC = get.partialbicoh(me)
           %Component of bicoherence attributed to the feature
          FF = me.wavefft(me.freqindx.Is);
          FF(:,me.order) = conj(FF(:,me.order));
          FF = prod(FF,2);
          FF(end+1) =0;
          pBC = FF./me.D;
          pBC = pBC(me.freqindx.remap);
          pBC(me.freqindx.PDconj) = conj(pBC(me.freqindx.PDconj));
          
        end
        function out = get.BIAS(me)        
          bias = sqrt(me.BIASnum./(me.D.^2+eps));
          out = bias(me.freqindx.remap);
        end
        function set.BIAS(me,in)
             in(isnan(in))=0;
             if min(size(in))==1
                me.BIASnum = in*me.D.^2; 
            else
                me.BIASnum = [in(me.freqindx.reduce)*me.D.^2;0];
             end
        end
        function out = get.fullmap(me)
           out = me.freqindx.remap; 
        end
        function set.lowpass(me,a)
           me.lowpassval = a./me.sampling_rate;
            me.update_frequency_indexing;
        end
        function set.shighpass(me,a)
           me.shighpassval = a./me.sampling_rate;
            me.update_frequency_indexing;
        end
        function set.slowpass(me,a)
           me.slowpassval = a./me.sampling_rate;
            me.update_frequency_indexing;
        end
        function set.glowpass(me,a)
           me.glowpassval = a./me.sampling_rate;
            me.update_frequency_indexing;
        end
        function set.xlowpass(me,a)
           me.xlowpassval = a./me.sampling_rate;
            me.update_frequency_indexing;
        end
        function set.highpass(me,a)
           me.highpassval = a./me.sampling_rate;
            me.update_frequency_indexing;
        end
        function set.xhighpass(me,a)
           me.xhighpassval = a./me.sampling_rate;
            me.update_frequency_indexing;
        end
        function out = get.lowpass(me)
          out = me.lowpassval*me.sampling_rate ;
        end
        function out = get.xlowpass(me)
           out = me.xlowpassval*me.sampling_rate;
        end
        function out = get.glowpass(me)
           out = me.glowpassval*me.sampling_rate;
        end
         function out = get.shighpass(me)
           out = me.shighpassval*me.sampling_rate;
         end
         function out = get.slowpass(me)
           out = me.slowpassval*me.sampling_rate;
        end
        function out = get.highpass(me)
          out =  me.highpassval*me.sampling_rate;
        end
        function out = get.xhighpass(me)
           out = me.xhighpassval*me.sampling_rate;
        end
        function out = get.buffersize(me)
            out = me.bufferN;
        end
        function set.buffersize(me,N)
           
           Norig = me.bufferN;
           me.bufferN = N;
           if islogical(me.padN)&& me.padN
               me.fftN = me.bufferN + N;
           elseif me.padN<1 && me.padN>0
               me.fftN = me.bufferN + round(me.padN*N);
           else
               me.fftN = me.padN + me.bufferN;
           end
           me.win = window(me.wintype,N); 
           me.radw = ifftshift((0:me.fftN - 1 )' - floor((me.fftN)/2))/(me.fftN)*2*pi;
           me.sampt = ifftshift((0:me.fftN - 1 ) - floor((me.fftN)/2)'); 
           if (N ~= Norig || length(me.freqs{1})~=me.fftN) && me.do_indexing_update
                freqs = {fftfreq(me.fftN)*me.sampling_rate};
                freqs(1:me.order) = freqs;
                me.freqs = freqs;
                me.update_frequency_indexing;
           end
           me.reset;
        end
        function out = get.window(me)
           out =  me.wintype;
        end
        function set.window(me,win)
           me.wintype=win;
            me.win = window(win,me.bufferN); %#ok<*CPROPLC>
        end
        function out = get.filterftlag(me)
           %%% Filter FT without circular shift adjustment
            out = zeros(me.fftN,size(me.G,2),size(me.G,3));
            out(me.keepfreqs{1},:,:) = me.G;
        end
        function out = get.filterfft(me)
           
            %%% Filter with lag adjustment
            out =me.filterftlag;
           
            %Adjust centering.
           dt = atan2(imag(me.lag),real(me.lag))/(2*pi)*me.fftN;
           delt = me.radw*dt;
           if isempty(me.radw)
               delt = 0;
           end
           out= exp(-1i*delt).*out;

        end
        function set.filterfft(me,in)
           
           in(isnan(in))=0;
           dt = atan2(imag(me.lag),real(me.lag))/(2*pi)*me.fftN;
           delt = me.radw*dt;
           F= exp(1i*delt).*in;
           
            me.G = F(me.keepfreqs{1},:,:) ;
            
        end
        function out = get.filterfun(me)
           %%% Filter function with lag adjustment
            F = me.filterfft;
%             out = ifftshift(real(ifft(F)));
            out = real(ifft(F));
            
        end
         function out = get.wavefft(me)
           
            F = me.waveftlag;
        %    [~,mxi] = max(ifft(F.*me.filterftlag));
            
               %Adjust centering.
            dt = atan2(imag(me.lag),real(me.lag))/(2*pi)*me.fftN;% + me.sampt(mxi);
            delt = me.radw*dt;
            out = exp(1i*delt).*F;
         
%             out = real(ifft(F));
            
        end
         function set.wavefft(me,in)
           
             F = in;
           
            F(isnan(F))=0;
            %Adjust centering.
            dt = atan2(imag(me.lag),real(me.lag))/(2*pi)*me.fftN;
            delt = me.radw*dt;
            if isempty(delt)
                delt = 0;
            end
            if ~isempty(F)
                me.waveftlag = exp(-1i*delt).*F;
            end        
        end
        function out = get.waveform(me)
           
            out = real(ifft(me.wavefft));
                                 
        end
        function set.waveform(me,in)
           
            if ~isempty(in)
               in(isnan(in))=0; 
               F = fft(in);
               me.wavefft = F;
            end   
        end
        
        function set.filterfun(me,in)
           
            if length(in)>me.bufferN
                warning('Filter function size does not match current buffer. Filter will be truncated.')
               in(me.bufferN+1:end)=[];
            elseif length(in)<me.bufferN
                warning('Filter function size does not match current buffer. Filter will be padded.')
                in(end+1:me.bufferN) = 0;
            end
            in(isnan(in))=0;
            in(end+1:me.fftN) = 0;
%             F =fft(fftshift(in));
            F =fft((in));
            me.filterfft = F;
            
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
%             out = ifftshift(me.waveform); 
             out = me.waveform;
        end
         %%%%%%%
        function set.feature(me,in)
%            me.waveform = fftshift(in); 
           me.waveform = in; 
        end
        %%%%%%%%
        function [Xfilt,FXshift,sgn] = apply_filter(me,X,apply_window,return_shifted,varargin)
            if nargin<3 || isempty(apply_window)
                apply_window = true;
            end
            if nargin < 4 || isempty(return_shifted)
               return_shifted = true; 
            end           
%             if nargin < 5 || isempty(center_delays)
%                center_delays = false; 
%             end
            FXshift = [];
            sgn = 1;
            if isempty(me.sampweight)
               smpw =1 ;
            else
                smpw = me.sampweight;
            end
            if size(X,1) == me.bufferN
                if apply_window
                    win = me.win;
                else
                    win =ones(size(X,1),1);
                end
%                 Xwin = fftshift(repmat(win,1,size(X,2)).*X,1);
                Xwin = repmat(win,1,size(X,2)).*X;
                Xwin(end+1:me.fftN,:) = 0;
                FXwin = fft(Xwin);
%                 FXwin = fft(X)';
                Xfilt = real(ifft(FXwin.*repmat(me.filterfft,1,size(X,2))));   
                if isscalar(smpw)
                    [~,mxi] = max(Xfilt.^me.order);
                else
                    [~,mxi] = max(Xfilt.^me.order.*repmat(smpw',size(Xfilt,1),1));
                end
                if mod(me.order,2)==0 
                   sgn = sign(Xfilt(mxi + (0:size(Xfilt,2)-1)*size(Xfilt,1)).*smpw'); 
                end
                if nargout >1 && return_shifted
                 
    %                 FX = fft(X);
                    samptc=(me.sampt);
                     dt = samptc(mxi);
        
                     me.delay = dt;             
%                      %%% Center the delays
%                      if center_delays
%                          
%                          mph = me.avg_delay;
%                          mdt = atan2(imag(mph),real(mph))*me.bufferN/(2*pi);
%                          dt = dt-mdt;
%                      end
                    delt = me.radw*dt;
                    FXshift = exp(1i*delt).*FXwin;
       
                else
                     FXshift = FXwin;
                end
                FXshift = FXshift.*repmat(sgn,size(FXshift,1),1);
%                 FXshift = FXshift*diag(sgn);
            else
%                  Xin = X(:);
                Xin = X;
                Xin(end+me.fftN,:) = 0;
                Xfilt = filter(me.filterfun,1,Xin);
                Xfilt = Xfilt(ceil(me.fftN/2)+1:end-floor(me.fftN/2),:);
                
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
            in(isnan(in))=0;
            if min(size(in))==1
                me.Bval = in; 
            else
                me.Bval = [in(me.freqindx.reduce);0];
            end
        end
        function set.D(me,in)
            in(isnan(in))=0;
            if min(size(in))==1
              	me.Dval=in; 
            else
                me.Dval  = [in(me.freqindx.reduce);0];
            end
        end
        function  set.Bpart(me,in)          
            if isnumeric(in)
                me.Bpartval={};
            
            elseif isempty(in) || min(size(in{1}))<=1
               	me.Bpartval=in; 
            else
                for kk = 1:length(in)
                    me.Bpartval{kk} = [in{kk}(me.freqindx.reduce);0];
                end
            end
        end
        %%%%%%%
        function out = get.H(me)
            
%             B = me.B(me.freqindx.remap); %#ok<*PROP>
%             B(me.freqindx.PDconj) = conj(B(me.freqindx.PDconj));
%             out = conj(B)./me.D(me.freqindx.remap).^2;
           BC = me.B./(me.D+eps);
           bias = sqrt(me.BIASnum./(me.D.^2+eps));
           bias(isnan(bias))=0;
           BC = (abs(BC)-bias).*BC./(abs(BC)+eps);
           H = BC./(me.D+eps);
           H = H(me.freqindx.remap);
           H(me.freqindx.PDconj) = conj(H(me.freqindx.PDconj));
           out = conj(H);
        end
        
        %%%%%%%
        function set.regressor(me,xin)
            if isa(xin,'regressor')
                if isempty(xin)
                    return
                end
                R = xin;
                xin = R.value;
            else
                R = [];
            end
            if size(xin,1)> me(1).buffersize
                for k = 1:size(xin,2)
                    X(:,k) = sum(me(1).chop_input(xin(:,k)),1);
                end
            elseif size(xin,1)==me(1).buffersize
                for k = 1:size(xin,3)
                    X(:,k) = sum(xin(:,:,k));
                end
            else
                X = xin;
            end
            if isempty(R)
                if isempty(which('regressor'))
                    return
                else
                    R = regressor(X);
                end
            else
                R2 = regressor(X);
                R.value = R2.value;
                if isfield(R,'m')
                    R.m=[];
                end
%                 R.m = R2.m;
            end
            for k = 1:length(me)
                me(k).regval = R;
            end
        end
        function out = get.regressor(me)
            out = me.regval;
        end
        
        function factreg(me,Fin,varargin)
            
            R = fact2reg(Fin,'ignore',Fin==0,varargin{:});
                      
            me.regressor = R;
        end
        %%%%%%%
        function out = xfilt(me,in,apply_window)
           if nargin < 2
               in = me.dat;
           end
           if nargin < 3 || isempty(apply_window)
              apply_window = false; 
           end
            [out,~] = me(1).apply_filter(in,apply_window,false);  
        
%             if size(in,1)==me(1).buffersize %% Make sure the output is consistent if the input happens to be of buffersize length
%                 out = ifftshift(out,1);
%             end
            
           if length(me)>1
               out = cat(sum(size(in)>1)+1,out,me(2:end).xfilt(in-me(1).xrec(in,[],apply_window),apply_window));
           end
        end
        %%%%%%%
        function out = ximp(me,in,return_sparse,apply_window)
           % Just get the thresholded times (for backward compatibility)
           if nargin < 2
               in = me.dat;
           end
           if nargin < 3  || isempty(return_sparse)
               %Cannot do sparse output if it requires more than 2
               %dimensions
               return_sparse = size(in,2) < 2 || length(me) < 2;
           end
           if nargin < 4 || isempty(apply_window)
              apply_window = false; 
           end
           out = me.xthresh(in,return_sparse,apply_window)>0;  
%            if length(me)>1
%                out = [out,me(2:end).ximp(in-me(1).xrec(in))];
%            end         

        end
        function out = xthresh(me,in,return_sparse,apply_window)
             % Get the thresholded data 
           if nargin < 2 || isempty(in)
               in = me.dat;
           end
           if nargin < 4 || isempty(apply_window)
              apply_window = false; 
           end
           if nargin < 3 || isempty(return_sparse)
               %Cannot do sparse output if it requires more than 2
               %dimensions
               return_sparse = size(in,2) < 2 || length(me) < 2;
           end
           
           xf = me(1).apply_filter(in,apply_window,false);
%            if size(in,1)==me(1).buffersize %% Make sure the output is consistent if the input happens to be of buffersize length
%                xf = ifftshift(xf,1);
%            end

           
           if return_sparse
               out = sparse(double(me(1).filter_threshold(xf)));  
           else         
               out = double(me(1).filter_threshold(xf));       
           end
           
            
           if length(me)>1
               out= cat(sum(size(in)>1)+1,out,me(2:end).xthresh(in-me(1).xrec(in,[],apply_window),return_sparse,apply_window));
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
           out(isnan(out)) = 0;
           if length(me)>1
               out =  cat(sum(size(in)>1)+1,out,me(2:end).xrec(in-out,thresh,apply_window,varargin{:}));
           elseif all(out(:)==0)
               sz = num2cell(size(in));
               sz{sum(size(in)>1)+1}=length(me);
               out(sz{:})=0;
           end
        end
        %%%%%%%
        function FFXpart = update_bispectrum(me,FXs,initialize)
            
            %Right now this updates in chunks with temporal decay weighting
            %applied only serially. That is, a simple average is obtained
            %for each chunk, which is then 
            
            if ~iscell(FXs)
                FXs  = {FXs};
            end
            
            m = size(FXs{1},2);
            
            if nargin < 3 || isempty(initialize)
                initialize = false;
            end
                
               
            if isempty(FXs{1})
                return
            end
            
               %%% Adjust for lag
            dt = atan2(imag(me.lag),real(me.lag))/(2*pi)*me.fftN;
            delt = me.radw*dt;
            delt(isnan(delt))=0;
            for k = 1:length(FXs)
                FXs{k} = repmat(exp(-1i*delt),1,size(FXs{k},2)).*FXs{k};
            end
            
            if length(FXs) == me.order
                FX = FXs{me.order};                
            elseif length(FXs)>1 && length(FXs)<me.order
                FXs = [FXs(ones(1,me.order-length(FXs)-1)),FXs,FXs(1)]; 
                FX = FXs{1};
            else
                FX = FXs{1};
            end
         
%             Xwin = Xin.*me.win;
%             FX = fft(Xwin);
           % FFX = 1;
%            FFXpart = ones([size(me.freqindx.Is,1),size(FX,2),me.order]);
            FFX = conj(FX(me.freqindx.Is(:,me.order),:));
            FFXpart = {};
            FFXpart(1:me.order-1) = {FFX};
            FFXpart{me.order} = ones(size(FFX));
            
            for k = me.order-1:-1:1
                if k>length(FXs)
                    FX = FXs{1};
                else
                    FX = FXs{k};
                end
                FXk = FX(me.freqindx.Is(:,k),:);
                for kk = setdiff(1:me.order,k) %%% Need multiple me.order symmetry regions for avg. partial cross-polyspectra
                    FFXpart{kk} = FFXpart{kk}.*FXk;
                end
                FFX = FFX.*FXk;
            end
            
            if isempty(me.sampweight)
                wgt = ones(size(FFX,2),1)/size(FFX,2);
            else
                wgt = me.sampweight;
            end
            
%             BX = mean(FFX,2);
            BX = FFX*wgt;
            XPSD = mean(abs(FX).^2,2);
%             BXpart = mean(FFXpart,2);
            
            BX(end+1,1) = 0;
            XPSD(end+1,:) =0;
%             BXpart(end+1,:) = 0;
            BXpart = {};
            for kk = 1:me.order
%                BXpart{kk} = mean(FFXpart{kk},2); 
               BXpart{kk} = FFXpart{kk}*wgt; 
               BXpart{kk}(end+1,:) = 0;
            end
            
            if initialize
                lradj = 1;
                fflr = 1;
                lrbias = 1;
                lr=1;
                me.sumlr = 1;
                me.sumlr2 = 1/size(FX,2);

            else
                [lradj,lr] = me.learningfunction(me.hos_learning_rate,m);
                 fflr = me.learningfunction(me.filter_adaptation_rate,m,1./me.filter_adaptation_rate);
%               fflr = (1-(1-me.filter_adaptation_rate)^m);
                %%% Adjust the learning rate according to the number of samples
                %%% in Xin, giving the sample weight the same total as if it
                %%% had been added serially.
                 asympedf = 2./lr-1; %Asymptotic EDF

                lrbias = 1./asympedf*(1-(1-lr).^(2*m)); % "Learning rate" for the sum of squared weights in the bias term
                me.sumlr = me.sumlr*(1-lradj) + lradj;
                me.sumlr2 = me.sumlr2*(1-lr).^(2*m) + lrbias;
            end            
            BX(isnan(BX))=0;
            me.B = (1-lradj)*me.B + lradj*BX;
%             me.Bpart = (1-fflr)*me.Bpart + fflr*BXpart;
             for kk = 1:me.order
                me.Bpart{kk} = (1-fflr)*me.Bpart{kk} + fflr*BXpart{kk};
             end
            XPSD(isnan(XPSD))=0;
            me.PSD = (1-lradj)*me.PSD + lradj*XPSD;
%            me.sumlr2 = (me.sumlr2-1./asympedf)*(1-me.current_learning_rate).^(2*m) + lrbias;
            
            switch me.normalization
                case 'awplv'
%                     NX = mean(abs(FFX),2);
                    NX = abs(FFX)*abs(wgt);
                    NX(end+1,1) = 0;
%                     XbiasNum = sum(abs(FFX).^2,2)./m^2;
                    XbiasNum = abs(FFX).^2*abs(wgt).^2;
                    XbiasNum(end+1,1) = 0;
                    me.BIASnum = me.BIASnum.*(1-lr).^(2*m) + lrbias*XbiasNum;
                    me.D = (1-lradj)*me.D + lradj*NX+eps;
                   
                case {'bicoh','bicoherence'}
                    %This doesn't seem to have strictly correct symmetry
                    XBCpart = mean(abs(FFXpart{1}).^2,2);
                    XBCpart(end+1,1) = 0;
                    me.BCpart = (1-lradj)*me.BCpart + lradj*XBCpart;                    
                    me.D = sqrt(me.BCpart.*me.PSD(me.freqindx.Is(:,1)))+eps;
            end
         
        end
        
         %%%%%%%
        function Gout = partial_delay_filt(me,Xs,returnfull,use_sample_bispectrum)
            
            
            
            if nargin < 4 || isempty(use_sample_bispectrum)
                use_sample_bispectrum = false; % Uses precomputed statistics if false
            end
            if nargin < 3 || isempty(returnfull)
                returnfull = true;
            end
            if ~iscell(Xs)
                Xs  = {Xs};
            end
            
            %m = size(Xs{1},2);
            
                               
            if isempty(Xs{1})
                return
            end
            
               %%% Adjust for lag
%             dt = atan2(imag(me.lag),real(me.lag))/(2*pi)*me.fftN;
%             delt = me.radw*dt;
%             delt(isnan(delt))=0;
%             delt = 0;
            for k = 1:length(Xs)
%                 FXs{k} = repmat(exp(-1i*delt),1,size(Xs{k},2)).*fft(Xs{k});
                FXs{k} = fft(Xs{k});
            end
            
            if length(FXs) == me.order
                FX = FXs{me.order};                
            elseif length(FXs)>1 && length(FXs)<me.order
                FXs = [FXs(ones(1,me.order-length(FXs)-1)),FXs,FXs(1)]; 
                FX = FXs{1};
            else
                FX = FXs{1};
            end
         
%             Xwin = Xin.*me.win;
%             FX = fft(Xwin);
           % FFX = 1;
%            FFXpart = ones([size(me.freqindx.Is,1),size(FX,2),me.order]);
            
            if ~use_sample_bispectrum
                FFX = conj(FX(me.freqindx.Is(:,me.order),:));
                FFXpart = {};
                FFXpart(1:me.order-1) = {FFX};
                FFXpart{me.order} = ones(size(FFX));

                for k = me.order-1:-1:1
                    if k>length(FXs)
                        FX = FXs{1};
                    else
                        FX = FXs{k};
                    end
    %                 FX(isnan(FX)) = 0;
                    FXk = FX(me.freqindx.Is(:,k),:);
                    for kk = setdiff(1:me.order,k) %%% Need multiple me.order symmetry regions for avg. partial cross-polyspectra
                        FFXpart{kk} = FFXpart{kk}.*FXk;
                    end
                    FFX = FFX.*FXk;
                end
            else
                FFXpart = me.update_bispectrum(FX,true);
            end
            
            BC = me.B./(me.D+eps);
            bias = sqrt(me.BIASnum./(me.D.^2+eps));
            bias(isnan(bias))=0;
            BC = (abs(BC)-bias).*BC./(abs(BC)+eps);
            H = conj(BC./(me.D+eps));
           
            Gpart = 0;
%             Bcheck = 0;
%              HFcheck =  H(1:end-1).*FFX; %Sanity check to make sure the integration is correct.
            len = length(me.B)-1;
            for k = 1:length(FFXpart)
                
                if k<=length(me.Imats)
                    I = me.Imats{k};
                    Iconj = me.Iconjmats{k};
                else
                    Iconj = integrator(me.freqindx.remap.*cast(me.freqindx.PDconj & me.freqindx.partialSymmetryRegions==k,class(me.freqindx.remap)),1,len,[0 len+1]);
                    I = integrator(me.freqindx.remap.*cast(~me.freqindx.PDconj & me.freqindx.partialSymmetryRegions==k,class(me.freqindx.remap)),1,len,[0 len+1]);
                    me.Imats{k}= I;
                    me.Iconjmats{k} = Iconj;
                end
                HF = repmat(H(1:end-1),1,size(FFXpart{k},2)).*FFXpart{k};
               
                Gpart = Gpart + I*HF + Iconj*conj(HF);
%                 Bcheck = Bcheck + I*HFcheck + Iconj*conj(HFcheck);
            end
            if returnfull
                Gout = zeros(size(Xs{1}));
                Gout(me.keepfreqs{1},:) = Gpart(me.keepfreqs{1}(abs(me.freqs{1})<=me.lowpass(1)),:);
            else
                Gout = Gpart(me.keepfreqs{1}(abs(me.freqs{1})<=me.lowpass(1)),:);
            end
        end
        
        %%%%%
        function update_filter(me)
    
%               postwin = @(x) 1+cos(2*pi*x);
%                postwin = @(x) exp(-x.^2*8);

%                Bpart = me.Bpart(me.freqindx.remap);

                Bpart = zeros(size(me.freqindx.remap));
                for k = 1:length(me.Bpart)                    
                    Bpart = Bpart + me.Bpart{k}(me.freqindx.remap).*(me.freqindx.partialSymmetryRegions==k);
                end
                Bpart(me.freqindx.PDconj) = conj(Bpart(me.freqindx.PDconj));

                GG = Bpart.*me.H;
            
                %%% Preserve the time winwowing
%                ng = size(GG);
%                t1 = (fftshift((0:ng(1)-1)-ceil(ng(1)/2)))./ng(1);
%                t2 = (fftshift((0:ng(2)-1)-ceil(ng(2)/2)))./ng(2);
%                [T1,T2]= ndgrid(t1,t2);
%                TMW=zeros(size(GG));
%                TMW(:) = postwin(T1).*postwin(T2).*postwin(T1-T2);
%                GGwin = fftn(real(ifftn(GG)).*TMW);
%                G = sum(GGwin,2);
               GG(isnan(GG))=0;
               G = sum(GG(:,:),2);

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
               
                if me.adjust_lag
                   ffun = ifftshift(real(ifft(me.filterftlag.*abs(me.waveftlag+eps))));                   
                   mph = sum(exp(-1i*2*pi*me.sampt(:)./me.fftN).*abs(ffun).^2)./sum(abs(ffun).^2);                   
                   mph = mph./(abs(mph)+eps);
                   if ~isnan(mph)
                      me.lag = mph; % Circularshift to keep filter energy centered on the window
                   end
               
%                    dt = atan2(imag(mph),real(mph))/(2*pi)*me.bufferN;
%                    delt = me.radw*dt;
%                     me.filterfft= exp(1i*delt).*me.filterfft;
                end
                
               %%% Also need to make sure that the output of the filter applied to
               %%% the feature waveform is centered with respect to the maximum!
               [~,mxi] = max(real(ifft(me.filterftlag.*me.waveftlag+eps)));
               if mxi~=1 && ~isnan(mxi)
                    me.filterfun = circshift(me.filterfun,-me.sampt(mxi));
               end
                   
        end
        
        function update_waveform(me,FXsh,initialize)
           m = size(FXsh,2);
           if nargin < 3 || isempty(initialize)
               initialize = false;
           end
           
           if initialize
               lradj=1;
           else
               lradj = me.learningfunction(me.filter_adaptation_rate,m,1./me.filter_adaptation_rate);
           end
               
%             lradj = (1-(1-me.filter_adaptation_rate)^m);

%            me.wavefft = me.wavefft*(1-lradj) + mean(FXsh,2)*lradj; 
           if isempty(me.sampweight)
               smpw = ones(size(FXsh,2),1)/size(FXsh,2);
           else
               smpw = me.sampweight;
           end
           me.wavefft = me.wavefft*(1-lradj) + (FXsh*smpw)*lradj; 
           
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
           % getsnip = min(me.bufferN-me.bufferPos,length(snip));
            me.inputbuffer(length(snip)+(1:end-length(snip))) = me.inputbuffer(1:me.bufferN-length(snip));
            me.inputbuffer(1:length(snip)) = snip;
            me.bufferPos = me.bufferPos + length(snip);
        end
        
        function out = get_input(me,xin,apply_window,use_shifted,initialize)
           
            % Process input if length is >= buffer size, else add to buffer.
            nxin = numel(xin);
              if nargin < 5 || isempty(initialize)
                 initialize=false;
              end
              if nargin < 4 || isempty(use_shifted)
                 use_shifted=true;
             end
             if nargin < 3 || isempty(apply_window)
                 apply_window=true;
             end
            if nxin >= me(1).bufferN
                me(1).bufferPos = 0; % Discard the buffer
                if  size(xin,1) ~=me(1).bufferN 
%                     stepn = round(me(1).poverlap*me(1).bufferN);
%                     nget = nxin - me(1).bufferN+1;
%                     tindx = (0:me(1).bufferN-1)';
%                     wint = (1:stepn:nget);
% 
%                     T = repmat(tindx,1,length(wint))+repmat(wint,length(tindx),1);
%                     Xchop = xin(T);
                    [Xchop,T] = me(1).chop_input(xin,false);
                    
                    snip = xin(T(end)+1:numel(xin));
                    if ~isempty(snip)
                        me(1).write_buffer(snip);
                    end
                else
                    Xchop = xin;
                end
                me(1).do_updates(Xchop,apply_window,use_shifted,initialize)
        
            else
                me(1).write_buffer(xin);
                return
            end    
            out = Xchop;
            if length(me)>1
               xrec = me(1).reconstruct(xin);
               xrec(isnan(xrec)&~isnan(xin))=0;
               
               out = [out,me(2:end).get_input(xin-xrec,apply_window,use_shifted,initialize)];
            end
            
            
            
        end
        function [Xchop,T] = chop_input(me,xin,apply_window,delay)
            if nargin < 3 || isempty(apply_window)
                apply_window = true;
            end
            if nargin < 4 || isempty(delay)
                delay=0;
            end
            nxin = length(xin);
            stepn = round(me(1).poverlap*me(1).bufferN);
            nget = nxin - me(1).bufferN+1;
            tindx = (0:me(1).bufferN-1)';
            wint = (1:stepn:nget)+delay;

            T = repmat(tindx,1,length(wint))+repmat(wint,length(tindx),1);
            T(T>length(xin))=length(xin);
            T(T<1)=length(xin);
            Xchop = xin(T);
            if apply_window
               Xchop = fftshift(repmat(me(1).win,1,size(T,2)).*Xchop,1);
            end
             
        end
        function  varargout = get_block(me,xin,maxiter,makeplot,segment,compno)
           
            % Fit a block of data all at once
            % Process input if length is >= buffer size, else add to buffer.
            if nargin < 6
                compno = 1;
            end
            Xsh=[];
            Xwin=[];
            T=[];
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
                if  size(xin,1) ~=me(1).bufferN 
                    if ~isfield(segment,'wint') || isempty(segment.wint)
                        stepn = round(me(1).poverlap*me(1).bufferN);
                        nget = nxin - me(1).bufferN+1;
                        wint = (1:stepn:nget);%./segment.fs;
                        segment.wint=wint/segment.fs;
                    else
                        wint = round(segment.wint*segment.fs);
                    end

                    T = repmat(tindx,1,length(wint))+repmat(wint,length(tindx),1);
                    T(T<1) = 1;
                    T(T>length(xin))=length(xin);
                    
%                     for k = 1:length(xin)
                        Xchop = xin(T);
                    
                        hasnans = any(xisnan(T));
                        segment.discarded = hasnans;
                        if all(hasnans)
                            fprintf('\nAll segments contain NaN values. Discarding these data')
                            if nargout > 1
                                 varargout = {Xsh,Xwin,T,wint,segment};
                            end
                            return
                        elseif any(hasnans)
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
                Xsh = Xwin;
                Xfilt=Xsh;
%                 end
                fprintf('\nComponent %3i Iter %3i',compno,0)
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
                std_moment = @(x)mean(cumulant(x,me(1).order,1)./(nanmean(x.^2).*nanmean(smpw.^2)).^(me(1).order/2));
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
                initialize = true;
                if me(1).use_partial_delay_method
                    me(1).get_input(Xsh,apply_window,use_shifted,initialize);
                    Gpart = me(1).partial_delay_filt(Xsh,false); 
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
                    elseif islogical(makeplot) && makeplot
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
                    fprintf('\b\b\b%03i',compno,k)
                    
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
                        me(1).feature = mean(Xsh,2);
                       end
                       if me(1).adjust_lag && me(1).do_filter_update
%                            ffun = ifftshift(real(ifft(me(1).filterftlag)));                   
                           ffun = ifftshift(real(ifft((me(1).filterftlag).*abs(me(1).waveftlag+eps))),1);                   
                           mph = sum(exp(-1i*2*pi*me(1).sampt(:)./me(1).fftN).*abs(ffun).^2)./sum(abs(ffun).^2);                   
                           mph = mph./(abs(mph)+eps);
                           if ~isnan(mph)
                             me(1).lag = mph; % Circularshift to keep filter energy centered on the window
                           end
                       end
                      
                   else
                       %%% This approach recomputes all statistics at every
                       %%% iteration, which is unnecessary.
                        me(1).get_input(Xsh,apply_window,use_shifted,initialize);
                        [Xfilt,FXsh] = me(1).apply_filter(Xsh,false,true);
                       newdt = me(1).delay;
                        Xsh = real(ifftshift(ifft(FXsh),1));
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
                %%% Also need to make sure that the output of the filter applied to
                %%% the feature waveform is centered with respect to the maximum!
                 [~,FXsh] = me(1).apply_filter(Xwin,false,true);
                 Xsh = real(ifft(FXsh));
                 me(1).feature = mean(Xsh,2);
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
                %%% Set the delays to the correct value for the original
                %%% data set;
               
                [~,~] = me(1).apply_filter(Xwin,apply_window);
                
                if length(segment.wint) == 0
%                     T = mod(repmat((0:size(Xsh,1)-1)',1,size(Xsh,2))+repmat(me(1).delay,size(Xsh,1),1),size(Xsh,1))+1;
%                     T = T+ones(size(T,1),1)*(0:size(T,2)-1)*size(T,1);
                else
                    T = T+ repmat(me(1).delay,size(T,1),1);
                    T(T<1)=1;
                    T(T>length(xin))=length(xin);

                    segment.wintadj = segment.wint+me(1).delay;
                end
            else
                me(1).write_buffer(xin);
                T=[];
            end    
            me(1).segment = segment;
            
            if length(me)>1
               xrec = me(1).reconstruct(xin);
               if nargout > 0
                   [Xsh2,Xwin2,T2] = me(2:end).get_block(xin-xrec,maxiter,makeplot,segment,compno+1);
                   Xsh = cat(3,Xsh,Xsh2);
                   Xwin = cat(3,Xwin,Xwin2);
                   T = cat(3,T,T2);
               else
                    me(2).Imats = me(1).Imats;
                    me(2).Iconjmats = me(1).Iconjmats;                    
                    me(2:end).get_block(xin-xrec,maxiter,makeplot,segment,compno+1);
               end
            end
            if nargout > 0
                varargout = {Xsh,Xwin,T,wint,segment};
            end
        end
        
        function out = hos_regress(me,yin,xin,varargin)
            
            % out = hos_regress(me,Y,X)
            % Linear regression on each element of the bicoherence matrix
            % computed for Y with regressors in X.
            %
            % see COMPLEXGLM
            
            do_permtest = false; %Do a permutation test to verify significant results
            maxpermn = 5e5; %#ok<NASGU>
            
            reg_args = {};
            if nargin < 3 || isempty(x)
                xin = ones(size(yin,1),1);
                reg_args = [reg_args,{'intercept',false}];
            end
            me(1).regressor = xin;
            
            if size(yin,2)==1
                Ychop = me(1).chop_input(yin,true);
                keep = ~any(isnan(Ychop)); %Discard any samples containing nans
                Ychop = Ychop(:,keep,:);
                for k = 1:size(xin,2)
                    Xchop = me(1).chop_input(xin(:,k),true);
                    Xchop = Xchop(:,keep,:);
                    x(:,k) = sum(Xchop)/sum(me(1).win);                
                end
                FY = fft(Ychop);
            else
                FY = fft(yin);
            end
            FFY = conj(FY(me(1).freqindx.Is(:,me(1).order),:));          
            for k = me(1).order-1:-1:1
                FYk = FY(me(1).freqindx.Is(:,k),:);
                FFY = FFY.*FYk;
            end
%             [~,out.dev0] = complexglm(FFX',ones(size(x,1),1),'diagonly',false,'intercept',false);
            [out.beta,out.dev,out.pval,out.iXX,out.sigma] = complexglm(FFY',x,'diagonly',false,reg_args{:});
            out.beta(:,end+1) = nan;
            out.dev(end+1) = nan;
            out.pval(end+1)=nan;
            
            out.sigma(end+1)=nan;
            me(1).regstat=out;
            se = sqrt(diag(out.iXX)*out.sigma);
            out.wald = out.beta./se;
            if do_permtest
               permi = 1;
               den = zeros(size(out.pval(1:end-1)));
               num = den;
               num_threshold = [5 4 3 2 1]; % Number of super-threshold permutations at which to stop at ceil(log10(permutation_index)
               fp = fprintf('\nPermutation %i, Nsig: %i',0,0);
               %%% Permutation test with number of permutations adjusted
               %%% to nominal p-value. 
               reseed;
               geti=1;
               while any(out.pval<=1/(permi*2)) && permi<maxpermn && sum(geti)>0
                   rp = randperm(size(FFY,2));
                   if mod(permi,10.^floor(log10(permi)-2))==0
	               geti = num < num_threshold(min(end,ceil(log10(permi+1))));
                       fp = fprintf([repmat('\b',1,fp),'Permutation %i, Nsig: %i'],permi,sum(geti))-fp;
                   end
                   [~,dev] = complexglm(FFY(geti,rp)',x,'diagonly',false);
                   den(geti) = den(geti)+1;
                   num(geti) =  num(geti)+(out.dev(geti)<dev);
                   permi=permi+1;
               end
               out.permp = num./den;
               out.permp(end+1) = nan;
%                out.nperms = den;
                
            end
         
%             if length(me)>1
%                out = [out,me(2:end).hos_regress(yin-me(1).xrec(yin),xin)];
%             end
         %             C = out.beta(:,1:end-1)*out.beta(:,1:end-1)';
%             [u,l] = svd(C);
%             me(1).regweight = u(:,1);
            
  %          out.tsq = sum((inv(out.iXX)*out.beta).*conj(out.beta))./out.sigma; 
        end
            
        
        
        function [Xthresh,Xcs,trialthresh] = filter_threshold(me,Xfilt,thresh,use_adaptive_threshold)
            
            % Apply a moment-based threshold
            if nargin < 3 || isempty(thresh)
                thresh= me.thresh;
            end
            if nargin < 4 || isempty(use_adaptive_threshold)
                use_adaptive_threshold = true;
            end
            
%             if size(Xfilt,2)==1
%                 Xcent = zscore(Xfilt(~isnan(Xfilt)));
%             else
                 zsc = @(x)(x-nanmean(x))./nanstd(x);
                 Xcent = zsc(Xfilt);
%             end
            %Xcent = Xfilt;
            if isempty(me.threshold_order)
                me.threshold_order = me.order;
            end
            Xmom = Xcent.^me.threshold_order;
            
            if size(Xfilt,1) == me.bufferN && size(Xfilt,2)==1 && use_adaptive_threshold
                 trialthresh = me.current_threshold;
                 Xcs = [];
            elseif me.threshold_order == 3
                % For the bispectrum compute normalized skewness
%                 keepsamples = ones(size(Xcent));
                  srt = sort(Xcent(:));
%                 outlier_threshold = 5;
                if me.outlier_threshold~=0
                   keepsamples = ~isnan(iterz(srt,me.outlier_threshold,-1)); % Suppress extreme negative outliers           
                else
                    keepsamples = ~isnan(srt);
                end
                srt(isnan(srt))=0;
                m1 = cumsum(srt.*keepsamples)./cumsum(keepsamples); % cumulative mean on sorted peaks
                m2 = cumsum(srt.^2.*keepsamples)./cumsum(keepsamples); % cumulative 2nd moment
                m3 = cumsum(srt.^3.*keepsamples)./cumsum(keepsamples); % cumulative 3rd moment
                %  Third cumulant
                c3 = m3 - 3*m2.*m1 + 2*m1.^3; % Third cumulant on sorted peaks
          
                keepsrt = srt>0 & c3>  thresh;
                detect = any(keepsrt);
                trialthresh = sum ((diff(keepsrt)>0).*srt(2:end,:)).^me.threshold_order;
                trialthresh(~detect) = Inf;
                Xcs=[];
            else
                % For now, apply a simple threshold on the standardized moment for
                % orders > 3. This should be improved to use the proper
                % cumulant.
                [Xsrt,srti] = sort(Xmom);
                Xpow = Xcent(srti).^2;
                keepsamples = ones(size(Xsrt));
                Mcs = cumsum(Xsrt.*keepsamples)./cumsum(keepsamples);
                Powcs = cumsum(Xpow.*keepsamples)./cumsum(keepsamples);
                if mod(me.threshold_order,2)==0
                    %%% Correction for power spectral component with even
                    %%% orders
                    Xbaseline = (me.threshold_order-1)*nanmean(Xcent.^2).^(me.threshold_order./2)+thresh;
                else
                    Xbaseline =thresh;
                end
                Mstd = Mcs./Powcs.^(me.threshold_order/2) - Xbaseline;
                Xthr = Mstd>thresh;
                Xthr = cumsum(diff([zeros(1,size(Xthr,2));Xthr])>0,'reverse')==0; % Use the last threshold crossing if there are multiple
                threshold_crossing = diff(Xthr)>0;
                detect =any(threshold_crossing);
                
                trialthresh(detect) = Xsrt(threshold_crossing)';
                trialthresh(~detect) = Inf;
            end
            
            if mod(me.threshold_order,2) < mod(me.order,2) %If using kurtosis for threshold setting with odd order HOSD, account for sign.
                Xmom(Xcent<0)=0;
            end
            
            switch me.threshold_type
                case 'hard'
                    THR = Xmom>=trialthresh;%repmat(trialthresh,size(Xmom,1),1);
                case 'soft'
                    sigm = @(x)1./(1+exp(-x));
                    THR = sigm(me.threshtemp*(Xmom-repmat(trialthresh,size(Xmom,1),1)));
                otherwise
                    error('Unrecognized threshold type')
            end
            Xthresh = Xfilt.*THR;
            Xthresh(isnan(Xthresh))=0;
        end
        
        %%%%%%
        function [Xrec,Xfilt] = reconstruct(me,X,threshold,apply_window)
            
            if nargin < 2
                X = me.inputbuffer;
            end
            if nargin < 3 || isempty(threshold)
                threshold = me.thresh;
            end
            if nargin < 4 
                apply_window = [];
            end
            xisnan = isnan(X);
%             X(xisnan)=0; 
            Xfilt = me.apply_filter(X,apply_window);
            if size(X,1) == me.bufferN
                 Xfilt = fftshift(Xfilt,1);
             
%                 win = me.win;
%               
%                 Xwin = repmat(win,1,size(X,2)).*X;

%             else
%                 Xwin = X;
            end

            % Xfilt = Xfilt(floor(me.bufferN/2):end-ceil(me.bufferN/2));
            Xthr=me.filter_threshold(Xfilt,threshold);
            
%             wf = fftshift(me.waveform);
            wf = me.waveform;
            wf(end+1:size(Xthr,1)) = 0;
            wf = circshift(wf,-floor(me.bufferN/2));
            Xthr(end+1:length(wf),:)=0;
            FXthresh =fft(Xthr);
            featfft = fft(wf);
            Xrec = real(ifft(FXthresh.*repmat(featfft,1,size(X,2))));
            Xrec(size(X,1)+1:length(wf),:) = [];
            %X(xisnan)=0;
            Xrec(xisnan)=0;

            use_filtered_lmse = true;
            if use_filtered_lmse
                %%% Apply the filter to the reconstructed data for LMSE fitting
                %%% so that the frequencies are appropriately weighted.
                Xrecfilt = me.xfilt(Xrec,apply_window);
                if size(X,1) == me.bufferN
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
            if nargin < 2
                me.reconbuffer = Xrec;
            end
        end
        function do_updates(me,X,apply_window,use_shifted,initialize)
            
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
            if me.do_filter_update
               me.update_filter; 
               Xsrt = mean(sort(zscore(Xfilt(:,getwin)).^me.order),2);
               lradj = me.learningfunction(me.filter_adaptation_rate,sum(getwin));
               me.thresholdbuffer = me.thresholdbuffer*(1-lradj) + lradj*cumsum(Xsrt);

            end
            if me.do_wave_update
              %  Xsh = real(ifft(FXsh(:,getwin)));
              
                me.update_waveform(FXsh(:,getwin),initialize); 
            end
            me.window_number = me.window_number+sum(getwin);
            me.outputbuffer = mean(Xfilt,2);
            me.shiftbuffer = real(ifft(mean(FXsh,2)));
           
        end
        
        function out = update_criteria(me,Xfilt) %#ok<INUSL>
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
        
        %%%%
%         function [T,Taligned] = chopmat(me,n)
%             
%             % Original and aligned, according to me.delay, chop matrices for a signal of length(n)
%             
%             if ~isscalar(n)
%                 n = length(n);
%             end
%             stepn = round(me(1).poverlap*me(1).bufferN);
%             nget = n - me(1).bufferN+1;
%             tindx = (0:me(1).bufferN-1)';
%             wint = (1:stepn:nget);
%             wintadj = wint + me.delay;
%             T = repmat(tindx,1,length(wint))+repmat(wint,length(tindx),1);
%             
%             
%         end
        function out = get.pad(me)
            out = me.padN;
        end
        function set.pad(me,in)
            me.padN = in;
            if islogical(me.padN) && me.padN
                me.fftN = me.bufferN+me.padN;
            else
                me.fftN = me.bufferN+me.padN;
            end
        end
    end
    
end

function out = fftfreq(N)
    out = ifftshift((0:N-1)-floor(N/2))/N;
end 
  
