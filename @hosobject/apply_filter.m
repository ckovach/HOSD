 function [Xfilt,FXshift,sgn] = apply_filter(me,X,apply_window,return_shifted,varargin)
        
% [Xfilt,FXshift,sgn] = apply_filter(me,X,apply_window,return_shifted)
%
% Apply the hosobject filter to a segment of data.
%
% Inputs:
%      X - data as a column vector or buffersize x N matrix
%      apply_window - Apply the window function specified in hos.window
%                     IF the input has length buffersize. Default is true.
%      return_shifted - Return the data shifted to the peak value in the
%                      second output argument if the length of X is buffersize.
%                       Default is true.
% Outputs:
%      Xfilt - Filtered data.
%      FXshift - FFT of filtered data, shifted to center at the peak value
%                (pos. peak for odd orders, |peak| for even) of Xfilt
%                if the input length is buffersize and return_shifted is 
%                true (otherwise it is unspecified). Also inverts according 
%                to the sign of the peak for even orders.
%      sgn     - Sign of the extremum for odd orders

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
    if size(X,1) == me.bufferN || size(X,1)==me.fftN
        if apply_window
            win = me.win;
            win(end+1:size(X,1))=0;
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