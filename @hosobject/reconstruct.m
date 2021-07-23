function [Xrec,Xfilt] = reconstruct(me,X,threshold,apply_window,use_adaptive_threshold)
      
% [Xrec,Xfilt] = reconstruct(me,X,[threshold],[apply_window],[use_adaptive_threshold])
%       
% Reconstruct a component signal according to the feature detection filter
% and waveform. This is done in the following steps:
%
%   1. Apply the feature detection filter to the input.
%   2. Threshold the filtered data such that the residual static cumulant of order
%       me.order is zero (or not greater than the value specified in
%       me.thresh).
%   3. Convolve the feature waveform with the output of step 2. 
%   4. Scale the result to minimize squared error. 
%
% Inputs:
%   X - input data as a column vector or me.buffersize x N array.
%   threshold - static cumulant threshold (default is me.thresh).
%   apply_window - If true and data has me.buffersize rows, applies window
%       specified in me.window. Default is set in
%       HOSOBJECT/APPLY_FILTER.
%   use_adaptive_threshold - If false, computes static cumulants from the
%       current sample, otherwise uses a pre-computed reference distribution.
%
% To obtain a complete multi-component decomposition use HOSOBJECT/XREC rather than
% HOSOBJECT/RECONSTRUCT.
%
% See also HOSOBJECT/CURRENT_THRESHOLD HOSOBJECT/FILTER_THRESHOLD
% HOSOBJECT/XREC HOSOBJECT/APPLY_FILTER
%
% Copyright Christopher K. Kovach, University of Iowa 2018-2021

    if nargin < 2
        X = me.inputbuffer;
    end
    if nargin < 3 || isempty(threshold)
        threshold = me.thresh;
    end
    if nargin < 4 
        apply_window = [];
    end
    if nargin < 5 || isempty(use_adaptive_threshold)
        use_adaptive_threshold=me.use_adaptive_threshold;
    end
    xisnan = isnan(X);
%             X(xisnan)=0; 
    Xfilt = me.apply_filter(X,apply_window);
    if size(X,1) == me.fftN %me.bufferN
         Xfilt = fftshift(Xfilt,1);

    end

    Xthr=me.filter_threshold(Xfilt,threshold,use_adaptive_threshold);

%             wf = fftshift(me.waveform);
    wf = me.waveform;
    wf(end+1:size(Xthr,1)) = 0;
    wf = circshift(wf,-floor(me.fftN/2));
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
        if size(X,1) == me.fftN
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