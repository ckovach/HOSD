function [out,snr] = xdetect(me,x); 

% xdet = xdetect(obj,x); 
%
% XDETECT detects instances of the feature in the input signal, x. It differs from
% XIMP and XTHRESH in that it applies a smoothing window whose width is determined from
% XTHRESH applied to the feature waveform. The reason for this is that some
% (in particular narrow band or periodic) features may generate multiple
% super-threshold peaks per feature instance, reflecting ambiguity in the
% timing of a a bandlimited or otherwise periodic feature. To avoid the inappropriate
% detection of multiple features instances in such caseses, XDETECT smooths
% the magnitude of the output of XTHRESH using a kernel of an appropriate scale for
% the feature and identifies features at the peak of the smoothed output.

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
    outnext = [];
    snrnext = [];
end


out = [pk==1,outnext];

if nargout > 1
    %Get the SNR of each peak based on root mean square power
    
    rmspow = sqrt(convn(abs(xthr).^2,g,'same'));

    snr = [rmspow.*(pk==1)/xfsd,snrnext];
end
%g = gausswin(12*smwin,6);


