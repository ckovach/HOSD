
classdef hosobject < handle
   
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
    % Copyright Christopher K. Kovach, University of Iowa 2018-2021
    
    properties
        
       % HOS order (3=bispectrum, 4=trispectrum, etc.)
       order = 3;
       
       % Axis frequencies
       freqs       
       
       % Power spectral density estimate
       PSD = 0;
      
       
       % Sampling rate of the input data
       sampling_rate = 1; 

       % Normalization used in computing bi-(or poly-)coherence
       normalization = 'awplv'; 
       
       %HOS learning rate for online mode
       hos_learning_rate = .01; 
       
       % Filter adaptation rate for online mode
       filter_adaptation_rate = .02; 
%        learningrate = .02; % Asymptotic learning rate
       filter_burnin = 20;
       hos_burnin = 20;
       
        % Nunmber of window processed
       window_number = 0; 
       
       % Default overlap between adjacent windows as the  proportion of step size to window width (poverlap = 1  is no overlap, poverlap = 2 interleaves a gap of 1 window duration between windows, etc.)
       poverlap = .5;    
      
       do_update = true;
       do_bsp_update = true;
       do_wave_update = true;
       do_filter_update = true;
       do_CDF_update = true;
       
       
       power_iterate = false; %Refine estimates with power iteration
       
       % Automatically apply a circular shift to the filter and waveforms to center the energy in both
       adjust_lag = true; 
       
       % A phasor representing the amount of circuglar shift added to the filter estimate according to angle (1 = no shift, +/-1i = max shift = buffersize/2)
       lag = 1; 
       thresh = 0;
       threshtemp = 1;
       threshold_type = 'hard';
       threshold_order = [];
       keepfreqs
       
       % Use only the principal domain in the estimates
       pdonly = true;
       dat = [];
       
       %For odd orders, reject negative outliers beyond this value. This improves the ability to recover features in the presence of noise and interference from similar features with opposite sign.
       outlier_threshold = 5;
       fftN = 1024;
       regstat = [];
       regweight=[];
       sampweight = [];
       
       mask = [];  %Mask applied to HOS domain. Only computes coefficients where the value is true;
       diagonal_slice = false; %Only compute coefficients on slices where at least one Wi==Wj, if true. If an integer, n, require at least n equalities.  
       %Matrices to integrate over all but 1st dimension of square form bicoherence, from the reduced form.
       Imats = {};
       Iconjmats = {};
       
       %Structure describing the segmentation of the input signal (if applicable)
       segment = struct('wint',[],'wintadj',[],'Trange',[],'fs',1,'discarded',[]);
       
       %(inverse) CDF is sampled at evenly spaced CDFupsample.*buffersize quantiles
       CDFupsample = 4; 
       
        %Buffer for the CDF of the filter output and moments up to order
       CDFbuffer=[];
        use_adaptive_threshold = true;
          running_mean=0;
        running_ssq = 1;
        running_var = 1;
        
        % During iteration use positive peaks only. Default true
        % if mod(order,2)=1 and false if mod(order,2) = 0.
        check_sign = false;
      
        %Include regions with these signatures, given as absolute
        %difference between the number of positive and negative
        %frequencies.
        include_signatures = []; 
        
        % Normalize each sample by its integrated magnitude spectrum to
        % suppress outliers.
        integrated_magnitude_normalization = false;
     end
  
    properties (GetAccess = public, SetAccess=protected)
       
        % Current input buffer
        inputbuffer = []; 
        
         % Current outputbuffer
        outputbuffer = [];
        reconbuffer = [];
        residualbuffer = [];
        shiftbuffer = [];
        
        % Structure containing information on frequency indexing
        freqindx = [];
        bufferPos = 0;
        sumlr=0;
        sumlr2=0;
        
        % Vector of sample frequencies in radian units
        radw = [];  
        
        % Vector of sample indices
        sampt = []; 
        
        % Vector of feature delays in the most recent input windows
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
       
        %%% If false, uses the older method which computed HOS  before
        %%% filer estimation, which is unnecessary in iterated estimation.
        %%% If true, then the partial_delay_filters are estimated for each
        %%% segment and the detection filter is obtained by shifting and 
        %%% averaging the segment filters. This option is included for the 
        %%% sake of debugging and will be removed eventually. 
         use_partial_delay_method = true;
       
       do_indexing_update = true;
    end
    properties (Dependent = true)
        
%         inputfft
%         inputfiltered  
%         input_shifted
        % Size of the input window
        buffersize
        
        % The feature waveform
        waveform; 
        %FFT of the feature waveform
        wavefft;
        
        % The ifftshifted filter function 
        filterfun %
        
        % FFT of the filterfunction
        filterfft %
        
        filterftlag;
        
        % Window applied to sample intervals
        window 
        
        %Bicoherence (or polycoherence) in square form
        bicoh 
        
        %Bicoherece (or polycoherence) in reduced form, as vector of non-redundant coefficients
        bicohreduced 
        partialbicoh
        H
        EDF
        highpass  ;
        
        % Lowpass on the edges (max freq); exclude regions in which the frequency axes exceeds this value
        lowpass ; %
        
        % Global lowpass: Exclude regions for which any frequency exceeds this value.
        glowpass ; %%
        
        % Highpass applied to sums of frequency pairs for order > 3.
        shighpass ; %%
        
        % Lowpass applied to sums of frequency pairs for order > 3.        
        slowpass ; %%
        
       % "OR" lowpass and highpass: include regions in which ANY of the frequencies meet the criterion 
       xlowpass ; 
       
       % "OR" highpass.This is useful to design filters selective for the interior or exterior regions of the bispectrum.
       xhighpass ;
       %Bispectrum (unnormalized) in square form
        Bfull 
        BIAS
        fullmap
        feature
        
        % Bispectrum (unnormalized) in reduced form, as vector of non-redundant coefficients.
        B 
        Bpart ;
        D ;
        pad; %
        regressor;
    end
    
    methods
       
       
        %%%% The following methods are defined within this file.
        function me = hosobject(order,varargin)
            % Hosobject contructor
            %
            %
            if nargin ==0
                return
            end
            if isa(order,mfilename) || isa(order,'hosminimal') || isa(order,'struct')
               obj = order;
%                fns = [{'BIASnum'};setdiff(properties(obj),{'BIAS','Bfull','H','bicoh','current_threshold','sampling_rate','freqindx','buffersize','filterftlag','fullmap','partialbicoh','filterfft','filterfun','bicohreduced'})];
              fns = setdiff(fieldnames(obj),{'freqs','BIAS','Bfull','H','bicoh','current_threshold','freqindx','filterftlag','fullmap','partialbicoh','bicohreduced'});
             
               metac = metaclass(me);
               props = metac.PropertyList;

               getprops = strcmp({props.SetAccess},'public');
               fns = [{'BIASnum'};intersect(fns,{props(getprops).Name})]; %This ensures that only fields with public set access are set to avoid unexpected behavior.
               
               if isa(order,'struct')
                   fnsunset = setdiff({'fftN','bufferN','sampling_rate','lowpass','freqs','freqindx'},fns);
                   for k = 1:length(fnsunset)
                       for kk = 1:length(obj)
                           obj(kk).(fnsunset{k})=me(1).(fnsunset{k});
                       end
                   end
               end               
               if length(obj)==1
                   obj(2:length(me)) = obj;
               elseif length(obj)>length(me)
                   me(length(obj)) = hosobject;
               end
               
               
               me(1).order = obj(1).order; 
               if obj(1).check_sign
                    me(1).check_sign = mod(obj(1).order,2)~=0;
               end
               me.initialize(obj(1).bufferN,obj(1).sampling_rate,obj(1).lowpass,obj(1).freqs,obj(1).freqindx,varargin{:})
               
               me(1).do_indexing_update = false;
               priority = intersect({'order','buffersize','pad','lag','sampling_rate','keepfreqs'},fns);% These fields should be set first
               for k = 1:length(priority)
                   me(1).(priority{k}) = obj(1).(priority{k});
               end
               fns = setdiff(fns,priority);
               for k = 1:length(fns)  
                   if isprop(me(1),fns{k})
                        try
                       me(1).(fns{k}) = obj(1).(fns{k});
%                        if me(1).lag~=obj(1).lag
%                            keyboard
%                        end
                        catch
                        end
                   end
               end
               me(1).do_indexing_update = true;
               if length(me)>1
                   me(2:end) = hosobject(obj(2:end));
               end
               return
            elseif  me(1).check_sign
               me(1).check_sign = mod(order,2)~=0;
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
            me(1).PSD = [z2;0];
%             me(1).G = ones(size(z));
            me(1).bufferPos = 0;
            me(1).B(:)=1;me(1).B=double(me(1).B);
            me(1).G(:)=1;me(1).G=double(me(1).G);
            me(1).D(:)=1; me(1).D=double(me(1).D);         
            me(1).BIASnum(:)=1; me(1).BIASnum = double(me(1).BIASnum);
            me(1).window_number=0;
            me(1).lag=1;
            me(1).Imats = {};
            me(1).Iconjmats = {};
            me(1).waveftlag=z;
            me(1).running_var = 1;
            me(1).running_ssq = 1;
            me(1).running_mean = 0;
            me(1).win = window(me(1).window,me(1).bufferN); %#ok<CPROP>
%             if me(1).fftN< me(1).bufferN || (~isempty(me(1).keepfreqs) && length(me(1).keepfreqs{1})~=me(1).bufferN)
%                 me(1).buffersize = me(1).bufferN;
% %                 me(1).fftN = me(1).bufferN;
%             end
            for k = 1:length(me(1).Bpart)
                me(1).Bpart{k}(:) = 0; me(1).Bpart{k} = double( me(1).Bpart{k});
            end
      
            me(1).waveform = z2;

            if isempty(me(1).threshold_order)
                me(1).threshold_order = me(1).order;
            end
            me(1).CDFbuffer=repmat(z2(:,ones(1,me(1).threshold_order)),me(1).CDFupsample,1);
            
            if length(me)>1
                me(2:end).reset();
            end
            
            %Set these to default values
            h = hosobject([]);
            fld =  {'do_bsp_update','do_bsp_update','do_filter_update','do_wave_update','do_update'};
            for k = 1:length(fld)
                me(1).(fld{k})=h.(fld{k});
            end            
        end
        
        
        function out = get.Bfull(me)
          out = me.B(me.freqindx.remap);
          out(me.freqindx.PDconj) = conj(out(me.freqindx.PDconj));
        end
        function BC = get.bicoh(me)
           %%% Bicoherence in square form
          BC = me.B./me.D;
           bias = sqrt(me.BIASnum./(me.D.^2+eps));
          BC = (abs(BC)-bias).*BC./(abs(BC)+eps);
          remap = me.freqindx.remap;
          BC(end+1:max(remap(:)),:)=nan;
      
          prm = circshift(1:me.order,-1);
          remap = remap + permute( cast(0:size(BC,me.order)-1,class(remap))',prm)*size(BC,1);
          BC = BC(remap);
          
          pdc = repmat(me.freqindx.PDconj,[ones(1,me.order-1) size(BC,me.order)]);
          BC(pdc) = conj(BC(pdc));
          
        end
        function BC = get.bicohreduced(me)
            %%% Bicoherence (or polycoherence) in reduced form (non-redundat coefficients as a
            %%% vector)
          BC = me.B./me.D;
           bias = sqrt(me.BIASnum./(me.D.^2+eps));
          BC = (abs(BC)-bias).*BC./(abs(BC)+eps);
          
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
        function pbcorr = pbcorr(me)
           %Correlation of feature partial polycohherence and sample polycoherence
          FF = me.wavefft(me.freqindx.Is);
          FF(:,me.order) = conj(FF(:,me.order));
          FF = prod(FF,2);
          FF(end+1) =0;
          pBC = FF./me.D;
          bc = me.bicohreduced;
          pbcorr = real((pBC'*bc)./sqrt(sum(abs(pBC).^2.).*sum(abs(bc).^2)));
        end
        function out = get.BIAS(me)        
          bias = sqrt(me.BIASnum./(me.D.^2+eps));
          out = bias(me.freqindx.remap);
        end
        function set.BIAS(me,in)
             in(isnan(in))=0;
             if min(size(in))==1
                me.BIASnum = in.*me.D.^2; 
            else
                me.BIASnum = [in(me.freqindx.reduce);0].*me.D.^2;
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
           if me.do_indexing_update && (N ~= Norig || length(me.freqs{1})~=me.fftN) 
                freqs = {me.fftfreq(me.fftN)*me.sampling_rate};
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
%             in(isnan(in))=0;
            if min(size(in))==1
                me.Bval = in; 
            else
                me.Bval = [in(me.freqindx.reduce);nan];
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
            setreg(me,xin);
        end
        function out = get.regressor(me)
            out = me.regval;
        end
        
        function factreg(me,Fin,varargin)
            
            R = fact2reg(Fin,'ignore',Fin==0,varargin{:});
                      
            me.regressor = R;
        end
        %%%%%%%
      
      
        %%%
       
        
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
      
        function out = update_criteria(me,Xfilt) %#ok<INUSL>
            out = true(1,size(Xfilt,2)); % Placeholder for now
        end
     
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
        
         %%%% These following lines declare methods that are defined in separate files within the
        %%%% @hosobject directory (see the files for explanations).
        update_frequency_indexing(me,freqindx,mask)
        %%%
        [lradj,lr] = learningfunction(me,learningrate,m,burnin)
        %%%
        [Xfilt,FXshift,sgn] = apply_filter(me,X,apply_window,return_shifted,varargin)
        %%%
        FFXpart = update_bispectrum(me,FXs,initialize)
        %%%        
        [Gout,sgn] = partial_delay_filt(me,Xs,returnfull,use_sample_bispectrum,normalization)
        %%%
        power_iteration(me,X)
        %%%
        varargout = get_block(me,xin,maxiter,makeplot,segment,compno,initialize)
        %%%
        out = get_input(me,xin,apply_window,use_shifted,initialize)
        %%%
        update_filter(me)
        %%%
        [Xchop,T,segment] = chop_input(me,xin,apply_window,delay,segment)
        %%%
        out = hos_regress(me,yin,xin,do_permtest,varargin)
        %%%
        [Xthresh,trialthresh] = filter_threshold(me,Xfilt,thresh,use_adaptive_threshold)
        %%%
        [Xrec,Xfilt] = reconstruct(me,X,threshold,apply_window,use_adaptive_threshold)
        %%%
        do_updates(me,X,apply_window,use_shifted,initialize)
        %%%
        out =current_threshold(me,Xcent,thresh)
        %%%
        out = xfilt(me,xin,apply_window)
        %%%
        out = xrec(me,xin,thresh,apply_window,varargin)
        %%
        out = xthresh(me,xin,return_sparse,apply_window)
        %%
        out = ximp(me,xin,return_sparse,apply_window)
        %%
        initialize(me,N,sampling_rate,lowpass,freqs,freqindex,varargin)

    end
    methods (Static)
        function out = fftfreq(N)
            out = ifftshift((0:N-1)-floor(N/2))/N;
        end 
    end    
    
end

  
