function [out,snr] = xdetect(me,x)

% [xdet,xsnr] = xdetect(obj,x); 
%
% XDETECT detects instances of the feature in the input signal, x. It differs from
% ximp and xthresh in that it applies a smoothing window whose width is determined from
% xthresh applied to the feature waveform  is that someestimate. The
% ratinale is that (in particular narrow band or periodic) features may generate multiple
% super-threshold peaks per feature instance, reflecting ambiguity in the
% timing of a single feature rather than multiple features. To address this, XDETECT smooths
% the magnitude of the output of XTHRESH using a kernel of an appropriate scale for
% the feature and assigns detections at the peaks of the smoothed output.
%
% Input:
%   obj - hosobject object.
%   x   - input signal.
% Output: 
%   xdet - logical array indicating whether a detection is present at each sample
%         of x.
%   xsnr - array giving the signal-to-noise ratio for each detaction according to RMS power of
%         xthresh within the smoothing window, normalized by the std dev. of
%         subthreshold samples of the filter output.
%
% See also XIMP, XTHRESH and XFILT

%C. Kovach 2025


% Get a characteristic duraion for the feature using fthresh
% This is a proxy for how the feature appears after filtering and
% thresholding in the original signal
fthr = me(1).xthresh(me(1).feature);

%Treat the magnitude of the thresholded feature like a probability
%distribution and compute the standard deviation.
wgt = abs(fthr);
wgt = wgt./sum(wgt);
smsd = sqrt(me(1).sampt.^2*wgt); %Standard deviation.

%Create a finite smoothing window with the same SD
g = hann(smsd.*5.5334);
g=g./sum(g);

%Get component, filter output and thresholded filter output
[xrec,xfilt,xthr] = me(1).xrec(x);


%Apply smoothing
xrsm = convn(abs(xthr),g,'same');


xfsd = nanstd(xfilt(~xthr));

[~,pk] = getpeak2(xrsm);

if length(me)>1
    if nargout == 1
         outnext = me(2:end).xdetect(x-xrec);
    else
        [outnext,snrnext] = me(2:end).xdetect(x-xrec);
    end
else
    outnext = logical([]);
    snrnext = [];
end


out = [pk==1,outnext];

if nargout > 1
    %Get the SNR of each peak based on root mean square power
    
    rmspow = sqrt(convn(abs(xthr).^2,g,'same'));

    snr = [rmspow.*(pk==1)/xfsd,snrnext];
end


