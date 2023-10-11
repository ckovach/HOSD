classdef ecgonline  < handle
    
    %Class for live feature identification, reconstruction and removal
    %using online HOSD. 
    %
    %To use:
    %
    %   1. Create object: ecg = ecgonline;
    %   2. Update with a data segment: ecg.update(x);
    % 
    % Input can be a column vector or a matrix, where the 2nd dimension is
    % channel. In the latter case, a multivariate version of HOSD will be
    % applied to estimate a spatio-temporal filter. 
    % 
    % Main properties:
    %   
    %   hos:     The hosobject implementing HOSD (see HOSOBJECT)
    %   feature: Current feature waveform estimate from HOSD.
    %   filterfun: Current detection filter estimate.
    %   xrec: Reconstructed component signal
    %   residual: Residual signal after removing xrec -- this is the
    %             "cleaned" signal.
    %   
    % See the commments in the script for additional properties. 
    %
    % See also HOSOBJECT
    
    %C. Kovach 2023
    
    properties
        hos      % hosobject
        input    % Input data segment
        xrec     % HOSD-based reconstructed component signal.
        xfilt    % Data filtered with estimated matched filter.
        filterfun  % Current matched filter estimate
        feature  % Current feature estimate
        residual % Residual after removing component (this is the "denoised" signal)
%        type = 'iterate'; %Apply the iterative algorithm - not working quite right yet
        type = 'stream'; %Use the streaming, non-iterative, algorithm.
        lowpass = 64; %Lowpass for hos estimation
        fs = 250; %Sampling rate
        hoswin = 8; %Analysis window duration in s
        buffersize = 7500; %Input duration in samples
        Nsamp = 0; %Cumulative number of segments sampled in HOSD estimation
        pre_highpass = 0.5 %Highpass before HOSD estimation
        pre_filter = []; % Filter function for highpass filter
        standardize = false; % Standardize each segment before estimation.
        learning_rate = 4e-3; %Learning rate for the bispectral running estimate
        outlier_threshold = 10; %Reject samples more than this number of standard deviations from the running mean
        running_average = []; % Running average, updated at the same rate as the HOSD filter
        running_mss = 1;    %Running mean square
        running_var = 1;    % Running variance (running_mms - running_average.^2)
        running_median = []; %Running average of sample medians;
        running_mad = 1; %Running average of absolute deviance from the median
    end
    
    methods
        
        function me = ecgonline(buffersize,fs,hosobj)
        
            if nargin < 1 || isempty(buffersize)
                me.buffersize = 7500;
            end
            if nargin < 2 || isempty(fs)
                me.fs = 250;
            end
            
            if nargin < 3 || isempty(hosobj)
                hosobj = hosobject(3,round(me.hoswin*me.fs),me.fs,me.lowpass);
            end
            hosobj.hos_learning_rate = me.learning_rate;
            hosobj.filter_adaptation_rate = me.learning_rate;
            hosobj.hos_burnin = 1./me.learning_rate;
%             hosobj.use_adaptive_threshold = true;
            %hosobj.poverlap = .75;
            me.hos = hosobj;
            
            me.pre_filter = fir1(2*floor(me.buffersize/10),me.pre_highpass/me.fs*2,'high');
            
        end
        
        function update(me,xin)
           
            %Pre-filter the data
            me.input = xin;

                        %Learning rate used by HOSD
            isn = any(isnan(xin),2);
         
            if isempty(me.running_average)
                me.running_average = nanmean(xin);
                me.running_median = nanmedian(xin);
            end
            if isnan( me.running_average), keyboard;end
     
            xin(isn,:) = repmat(me.running_average,sum(isn),1);
            xprefilt = filtfilt(me.pre_filter,1,xin);
            xprefilt(isn,:) = nan;
            
            x_for_thresholding = xin;
            
            if me.hos.EDF>10 
               %Reject outliers exceeding the specified threshold and
               %replace with nans
          
%                z = (xin-me.running_average)./sqrt(me.running_var);
               z = (x_for_thresholding-me.running_median)./me.running_mad;
               reject = any(abs(z)>me.outlier_threshold,2);
%                xprefilt(reject) = nan;
            %   if any(reject), keyboard;end
            else
                reject = false(size(xin,1),1); 
            end
            isn = isn | reject;
            if all(isn)
                return
            end
             if ~all(isnan(xin))
                lradj = me.hos.learningfunction(me.hos.filter_adaptation_rate,floor(sum(~any(isnan(x_for_thresholding),2))./me.hos.buffersize),me.hos.filter_burnin);            
                me.running_average = me.running_average*(1-lradj) + nanmean(x_for_thresholding + 0./~reject)*lradj;
                me.running_mss = me.running_mss*(1-lradj) + nanmean(x_for_thresholding.^2+ 0./~reject)*lradj;
                me.running_var = me.running_mss-me.running_average.^2;
      
                me.running_median = me.running_median*(1-lradj) + nanmedian(x_for_thresholding)*lradj;
                me.running_mad = me.running_mad*(1-lradj) + nanmedian(abs(x_for_thresholding-me.running_median))*lradj;
             end
           
            dx = xin-xprefilt; %The lowpass component will be added back in at the end
            xin = xprefilt;
            
            xin(isn,:) = nan;
           
            if me.standardize
                xsd = nanstd(xin);
                xm = nanmean(xin);
            else
                xsd = 1;
                xm = 0;
            end
            xin = (xin-xm)./xsd;
            x = xin + 0./~reject; %Replaces rejected values with nans here
    
            
            switch me.type
                case 'iterative'
                    me.hos.get_block(x,25,false,[],1,false);
                case 'stream'
                    me.hos.get_input(x);
            end
            
            xin(any(isn,2),:) = repmat(nanmean(xin),sum(any(isn,2)),1); %For the purpose of reconstruction, replace nans with the mean so we're less prone to miss events in the vicinity of nans.
            me.xrec = me.hos.xrec(xin)*xsd;
            me.residual = (xin-me.xrec)*xsd + xm + dx;
            me.xfilt = me.hos.xfilt(xin);
            me.feature = me.hos.feature;
            me.filterfun = me.hos.filterfun;
            me.Nsamp = me.hos.EDF;

        end
        function reset(me)
            me.Nsamp = 0;
            me.hos.reset();
        end
                
    end
end
