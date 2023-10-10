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
        hoswin = 2; %Analysis window duration in s
        buffersize = 7500; %Input duration in samples
        Nsamp = 0; %Cumulative number of segments sampled in HOSD estimation
        pre_highpass = 0.5 %Highpass before HOSD estimation
        pre_filter = []; % Filter function for highpass filter
        standardize = false; % Standardize each segment before estimation.
        learning_rate = 1e-3; %Learning rate for the bispectral running estimate
        outlier_threshold = 6; %Reject samples more than this number of standard deviations from the running mean
        running_average = 0; % Running average, updated at the same rate as the HOSD filter
        running_mss = 1;    %Running mean square
        running_var = 1;    % Running variance (running_mms - running_average.^2)
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
            hosobj.hos_burnin = 100;
%             hosobj.use_adaptive_threshold = true;
            %hosobj.poverlap = .75;
            me.hos = hosobj;
            
            me.pre_filter = fir1(2*floor(me.buffersize/10),me.pre_highpass/me.fs*2,'high');
            
        end
        
        function update(me,xin)
           
            %Pre-filter the data
            me.input = xin;

                        %Learning rate used by HOSD
            if ~all(isnan(xin))
                lradj = me.hos.learningfunction(me.hos.filter_adaptation_rate,floor(sum(~any(isnan(xin),2))./me.hos.buffersize),me.hos.filter_burnin);            
                me.running_average = me.running_average*(1-lradj) + nanmean(xin)*lradj;
                me.running_mss = me.running_mss*(1-lradj) + nanmean(xin.^2)*lradj;
                me.running_var = me.running_mss-me.running_average.^2;
            end
            isn = any(isnan(xin),2);
            xin(isn,:) = repmat(me.running_average,sum(isn),1);
            xprefilt = filtfilt(me.pre_filter,1,xin);
            
            
            if me.hos.EDF>10 
               %Reject outliers exceeding the specified threshold and
               %replace with nans
          
               z = (xin-me.running_average)./sqrt(me.running_var);
               reject = abs(z)>me.outlier_threshold;
%                xprefilt(reject) = nan;
            %   if any(reject), keyboard;end
            else
                reject = false; 
            end
            
            dx = xin-xprefilt; %The lowpass component will be added back in at the end
            xin = xprefilt;
            xin(isn) = nan;
           
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
