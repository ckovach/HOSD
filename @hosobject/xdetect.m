function [out,snr,xrsm,xfsd,xfilt,xthr,g] = xdetect(me,x,pow)

% [xdet,xsnr,xfsm,xfsd,xfilt,xthr,g] = xdetect(obj,x,pow); 
%
% XDETECT detects instances of the feature in the input signal, x. It differs from
% ximp and xthresh in that it applies a smoothing window whose width is determined from
% xthresh applied to the feature waveform, serving as a proxy for the timing ambiguity
% of individual feature instances. Detections are identified at peaks in the smoothed output.
% The rationale is that (in particular narrow band or periodic) features may generate multiple
% super-threshold peaks per feature instance, reflecting ambiguity in the
% timing of a single feature rather than multiple features. 
%
% Input:
%   obj - hosobject object.
%   x   - input signal.
%   pow - Compute SNR using an Lp norm of the given power. A higher value
%         returns an SNR that is closer the maximum amplitude of the
%         feature. Default is 10 (to approximate a max filter. 
% Output: 
%   xdet - logical array indicating whether a detection is present at each sample
%         of x.
%   xsnr - array giving the signal-to-noise ratio for each detaction according to RMS power of
%         xthresh within the smoothing window, normalized by the std dev. of
%         subthreshold samples of the filter output.
%   xrsm - smoothed thresholded output from whose peaks xdet is obtained.
%   xfsd - Estimated noise standard deviation in the detection filter output.
%   xfilt - Detection filter output.
%   xthr - Thresholded filter output.
%   g    - smoothing window applied to the xthr.

% See also XIMP, XTHRESH and XFILT

%C. Kovach 2025


if nargin < 3 || isempty(pow)
    pow = 2;
   % pow = 10; %Using root-mean-high-order-power so that the envelope is closer to the peak value as an approximation to the max filter, while keeping the output smooth
          %Hilbert amplitude smears too much in time and introduces
          %artifacts.
end

%Get component, filter output and thresholded filter output
[xrec,xfilt,xthr] = me(1).xrec(x);

% Get a characteristic duration for the feature using fthresh
% This is a proxy for how the feature appears after filtering and
% thresholding in the original signal

ff = me(1).xfilt(me(1).feature);
fthr = me(1).xthresh(me(1).feature);

%Treat the magnitude of the thresholded feature like a probability
%distribution and compute the standard deviation.
wgt = abs(fthr);
%%%Additionally, limit to the region with 99% of the total energy
toteng = ifftshift(cumsum(fftshift(abs(ff).^2)))./sum(abs(ff).^2);
wgt = wgt.*(toteng>.005 & toteng < .995);

%Normalize
wgt = wgt./sum(wgt);

smsd = sqrt(me(1).sampt.^2*wgt); %Standard deviation.

%Create a finite smoothing window with the same SD
nsamp= max(3,round(smsd.*pi./sqrt(pi^2/12-1/2)));
nsamp = nsamp + 1-mod(nsamp,2); %Keep the window size odd
g = hann(nsamp);
g=g./sum(g);


%xfsd = nanstd(xfilt(~xthr));

% xfsd = sqrt(nanmean(convn(xfilt.^2 + 0./~xthr,g,'same')));
xfsd = std(xfilt(~xthr));


%%% Adjust the sd estimate because it is estimated from a truncated
%%% distribution...
% thr = nthroot(me.current_threshold(xfilt./xfsd),me.order);
% if mod(me.order,2)==0;
%     lb = -thr;
% else
%     lb = -Inf;
% end
% ub = thr;
% sdadj = trnormcum(lb,ub,2);
% xfsd = xfsd./sqrt(sdadj);

%Apply smoothing
%xrsm = sqrt(convn(abs(hilbert(xthr)).^2,g,'same'));
%xrsm = nthroot(convn(abs(xthr./xfsd).^pow,g,'same'),pow);

%RMS normalized to estimated noise SD
xrsm = nthroot(convn((xthr./xfsd).^pow,g,'same'),2);
[~,pk] = getpeak2(xrsm);
%%% Peak weighted xrsms.
%%% Stabilize the denominator with 1 * max(g) so regularization scales with the kernel. 
norm = convn((xthr>0),g,'same')+1*max(g(:));
xmax = nthroot(xrsm.^pow./norm,pow);
%xmax = xrsm;


if length(me)>1
    if nargout == 1
         outnext = me(2:end).xdetect(x-squeeze(xrec),pow);
    else
        [outnext,snrnext,xrsmnext,xfsdnext,xfiltnext,xthrnext,gnext] = me(2:end).xdetect(x-squeeze(xrec),pow);
    end
else
    outnext = logical([]);
    snrnext = [];
    xrsmnext = [];
    xfsdnext = [];
    xfiltnext = [];
    xthrnext = [];
    gnext = [];
end


out = [pk==1 &xrsm>0,outnext];

if nargout > 1
    %Get the SNR of each peak based on root mean square power
    
    % rmspow = xrsm;
    rmspow = xmax; %max filtered data instead of rms to calculate peak snr
%    rmspow = sqrt(convn(abs(xthr).^2,g,'same'));

    snr = [rmspow.*(pk==1),snrnext];

    xrsm = [xrsm,xrsmnext];
    xfsd = [xfsd,xfsdnext];
    xfilt = [xfilt,xfiltnext];
    xthr = [xthr,xthrnext];
    if ~isempty(gnext)
        if size(gnext,1)< size(g,1)
            ln = size(gnext,1);
            gnext(size(g,1),:) = 0;
            gnext = fftshift(circshift(gnext,-floor(ln/2)),1);
        elseif size(gnext,1)>size(g,1)
            ln = size(g,1);
            g(size(gnext,1),:) = 0;
            g = fftshift(circshift(g,-floor(ln/2)),1);
        end
    end

    g = [g,gnext];
end


