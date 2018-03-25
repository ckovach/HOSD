

function [out,Xadj,X,dt]=bsident(x,segment,lpfilt,ncomp,varargin)

% out = bsident(x,segment,lpfilt,ncomp,varargin)
%
% Signal identification based on bispectral statistics. 
% The signal is first divided into overlapping windows and the delay
% correction implemented in BSTD is applied. The average of the
% delay-corrected result is used to identify the salient feature. 
%
% Input arguments:
%
%  x - the input signal.
%
%  segment - the duration of segment windows in samples. Alternatively can
%            be struct with details of more complicated segmentation schemes. 
%
%   lpfilt - lowpass filter as a proportion of the sampling freq. 
%
% Output arguments:
%
%   out - struct with the following fields:
%         .B     - Bicoherence
%         .wb    - Frequencies corresponding to rows and columns of B
%         .BFILT - Optimal filter for feature identification
%         .f     - feature identified in the average of the shift-corrected
%                  segments.
%         .xfilt - Signal filtered by BFILT
%         .xrec  - Reconstructed "denoised" signal obtained by thresholding
%                  xfilt (according to a simple threshold or kmeans) and 
%                  convolving with the feature,f. 
%         .dt    - delay correction for each segment following each iteration. 
%                  The final correction is  dt = sum(out.dt);
%         .ximp  - Indices of samples in xfilt that survive thresholding.
%         .segment - struct with the  details of the sementation:
%                      .Trange: time range of the segment window in samples
%                      .fs: Assumpled sampling rate (default = 1).
%                      .polvp: percent overlap of adjacent windows.
%                      .wint: location of segment windows (time of center).
%                           Segmentation matrix can be obtained with
%                           T = chopper(Trange,wint,1);
%                      .tt: window-relative time, i.e. Trange(1):1/fs:Trange(2)
%                      .wintadj: delay-adjusted window times.
%
% See also BSTD, CHOPPER
%
%
% C. Kovach 2017

niter = 10;

%%% Method to identify impulse
%impulse_method = 'zthreshold';
impulse_method = 'kmeans';
%decomp_method = 'residual';
decomp_method = 'direct';

type = 'mean';
% type = 'svd';

if nargin < 3 || isempty(lpfilt)
    lpfilt = .1;
end
if nargin <4 || isempty(ncomp)
    ncomp = 1;
end
% if nargin < 5 
%     regressors = [];
% end


n = length(x);
default_povlp=.5;
if isnumeric(segment)
    if isscalar(segment)
        segment= [-1 1]*segment/2;
    end
    segment = struct('Trange',segment,'fs',1,'povlp',default_povlp);
end    
if ~isfield(segment,'window')
    segment.window = @hann;
elseif isnumeric(segment.window)
    segment.window = @(varargin)segment.window;
end
    
if ~isfield(segment,'fs')
    segment.fs=1;
end

if ~isfield(segment,'wint') && isfield(segment,'povlp')
    segment.wint= 1/segment.fs:diff(segment.Trange)*(1-segment.povlp):n/segment.fs;
end

[T,tt] = chopper(segment.Trange, segment.wint,segment.fs);
segment.wint(:,any(T<1 | T>n)) = [];
T(:,any(T<1 | T>n)) = [];

segment.tt = tt;

nX = size(T,1);

out = struct('BFILT',[],'f',[],'dt',[],'xrec',[],'xfilt',[],'ximp',[],'a',[],'exvar',[],'B',[],'wb',[],'segment',[]);
xresid = x;
X = x(T);
thresh = 1;

for kk = 1:ncomp
%%
    switch decomp_method
        case 'residual'
            Fremove = [];
        case 'direct'
%            xresid=x;
            Fremove=[out.f];
    end
      X = xresid(T);
  

    Xadj = diag(segment.window(nX))*X;
    clear dt
    for k= 1:niter  %%% Apply bstd iteratively
        [dt(k,:),Xadj,B,BFILT,wb] = bstd(Xadj,lpfilt,Fremove,varargin{:});
        k
    end
    
    switch type
        case 'svd'
            [u,l] = svd(Xadj);

            [~,pckeep] = min(diff(diag(l)));
            f = u(:,1:pckeep)*sqrt(l(1:pckeep,1:pckeep));
        case 'mean'
            
            f = mean(Xadj,2);
    end
    
    f(n,:) = 0;
    f = circshift(f,-floor(nX/2));

    H = fftshift(BFILT);
    %H = H./sqrt(sum(abs(H).^2));
    H(nX) = 0;
    H = circshift(H,-floor(length(BFILT)/2));
    h = real(ifft(H));
    bfilt = h;
    h(n) = 0;
    h = circshift(h,-ceil(nX/2));

    xfilt = ifft(fft(xresid).*fft(h));

%     pk = getpeak(xfilt);
%     pk = pk(zscore(xfilt(pk))>thresh);
% 
%     ximp = zeros(size(x));
%     ximp(pk) = 1;

    %xrec = ifft(fft(ximp).*fft(f));%*nX/n;
    switch impulse_method
        case 'zthreshold'
            ximp = zscore(xfilt)>thresh;
        case 'kmeans'
            [km,kmc]  = kmeans(xfilt + 0./(zscore(xfilt)>1),2);
            [~,mxi] = max(kmc);
            ximp = km==mxi;
    end
    xrec = ifft(fft(xfilt.*ximp).*fft(f));%*nX/n;
    a = xrec'*xresid./sum(xrec.^2);
    xrec = a*xrec;
   
    out(kk).BFILT = bfilt;
    out(kk).f= mean(Xadj,2);
    out(kk).dt= sum(dt);
    out(kk).xrec = xrec;
    out(kk).xfilt = xfilt;
    out(kk).ximp= find(ximp);
    
    out(kk).a = a;
    out(kk).exvar = 1-sum(abs(xresid-xrec).^2)./sum(abs(xresid).^2);
    out(kk).B = B;
    out(kk).wb = wb;
    segment.wintadj = segment.wint+round(sum(dt))./segment.fs;
    out(kk).segment = segment;

     xresid =xresid-xrec;

end
    
    
