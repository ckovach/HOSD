classdef ecgonline  < handle
    
    %Class for online feature identification, reconstruction and removal
    properties
        hos      % hosobject
        input    % Input data segment
        xrec     % HOSD-based reconstructed component signal.
        xfilt    % Data filtered with estimated matched filter.
        filtfun  % Current matched filter estimate
        feature  % Current feature estimate
        residual % Residual after removing component (this is the "denoised" signal)
%        type = 'iterate'; %Apply the iterative algorithm - not working quite right yet
        type = 'stream'; %Use the streaming, non-iterative, algorithm.
        lowpass = 64; %Lowpass for hos estimation
        fs = 250; %Sampling rate
        hoswin = 2; %Analysis window duration in s
        buffersize = 7500; %Input duration in samples
%         M = 0; 
%         S = 0;
        Nsamp = 0; %Cumulative number of segments sampled in HOSD estimation
%         mahalthresh = 3; 
        pre_highpass = 0.5 %Highpass before HOSD estimation
        pre_filter = []; % Filter function for highpass filter
        standardize = false; % Standardize each segment before estimation.
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
            hosobj.hos_learning_rate = 1e-3;
            hosobj.filter_adaptation_rate = 1e-3;
            hosobj.hos_burnin = 1;
            %hosobj.poverlap = .75;
            me.hos = hosobj;
            
            me.pre_filter = fir1(2*floor(me.buffersize/10)+1,me.pre_highpass/me.fs*2,'high');
            
        end
        
        function update(me,xin)
           
            %Pre-filter the data
            xprefilt = filtfilt(me.pre_filter,1,xin);
           % dx = xin-xprefilt;
            xin = xprefilt;
             
            if me.standardize
                xsd = nanstd(xin);
                xm = nanmean(xin);
            else
                xsd = 1;
                xm = 0;
            end
            xin = (xin-xm)./xsd;
            x = xin;
            %x = me.hos.chop_input(xin,false);
            
            
            
%             nsamp = me.Nsamp+size(x,2);
%             m = me.M*me.Nsamp./nsamp + sum(x,2)./nsamp;
%             ss = me.S*me.Nsamp./nsamp + x*x'./nsamp;
%             
%             d = (sum((pinv(ss)*(x-m)).*(x-m),1));
%             
%             p = chi2cdf(d,size(x,1));
%             
%          %   keep = d<me.mahalthresh | nsamp < 2*size(x,1); %Reject outliers here
%             keep = p < .99;
%             
%             if any(~keep), keyboard, end
%             x = x(:,keep);
            nsamp = me.Nsamp + size(x,2);
            %Update mean and covariance after rejecting outliers
%             me.M = me.M*me.Nsamp./nsamp + sum(x,2)./nsamp;
%             me.S = me.S*me.Nsamp./nsamp + x*x'./nsamp;
            me.Nsamp = nsamp;
            
            
            switch me.type
                case 'iterative'
                    me.hos.get_block(x,25,false,[],1,false);
                case 'stream'
                    me.hos.get_input(x);
            end
            
            me.xrec = me.hos.xrec(xin)*xsd;
            me.residual = (xin-me.xrec)*xsd + xm;
            me.input = xin*xsd+xm;
            me.xfilt = me.hos.xfilt(xin);
            me.feature = me.hos.feature;
        end
        function reset(me)
            me.Nsamp = 0;
            me.hos.reset();
        end
                
    end
end
