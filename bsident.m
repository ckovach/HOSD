

function [out,Xadj,X,dt]=bsident(x,segment,lpfilt,ncomp,opts,varargin)

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

default_opts = struct('niter',10,...
'impulse_method','skew0',...%%% Method to identify impulses; 'skew0' retains samples such that remaining samples have 0 skewness
'impulse_skewness_sd_threshold',0,...%%% For the 'skew0' method, this sets the threshold in units of std dev. of the estimator.
'outlier_threshold',10,... %Exclude outliers above this value using ITERZ (iterated z-score) from the skewness computation
'decomp_method','residual',...%%% decompositions method
'resegment',false,... %If true, the signal is resegmented after each iteration. This allows segments to drift to any position within the signal.
'showprog',true,... %% Show a real-time plot of the realignment 
'use_ideal_filter',false,...
'pre_filter',true,... % Filter the signal so that it is zero-mean at the scale of the observation window
'skewness_threshold',-Inf,... %Exclude windows with filtered skewness below this threshold, under the assumption that do not contain the feature.
'lpfilt',.25,...
'ncomp',1); %Exclude windows with negative skewness after BFILT, assuming they do not contain the transient with high signal to noise ratio. 

if nargin > 4 && ~isempty(opts)
    if isstruct(opts)
        valid_fields = fieldnames(default_opts);
        fns = fieldnames(opts);
        for k = 1:length(fns)          
            if ismember(fns{k},valid_fields)
                default_opts.(fns{k}) = opts.(fns{k});
            else
                error('%s is not a recognized option.',fns{k});
            end
        end
    else
        varargin = {opts,varargin{:}};
    end
end
opts = default_opts;
% type = 'mean';
% type = 'svd';

if nargin > 2 &&  ~isempty(lpfilt)
    opts.lpfilt = lpfilt;
end
if nargin >3 && ~isempty(ncomp)
    opts.ncomp = ncomp;
end
% if nargin < 5 
%     regressors = [];
% end


n = length(x);
default_povlp=.5;
if isnumeric(segment)
    if isempty(segment)
        segment = opts.segment;
    elseif isscalar(segment)
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


wintorig = segment.wint;
[T0,tt] = chopper(segment.Trange, segment.wint,segment.fs);
segment.wint(:,any(T0<1 | T0>n)) = [];
T0(:,any(T0<1 | T0>n)) = [];

if ~isfield(opts,'segment')
        opts.segment = segment;
end

if ~isfield(opts,'bstdargs')
        opts.bstdargs = varargin;
end

segment.tt = tt;
nX = size(T0,1);


if opts.pre_filter
   
    b = fir1(nX,4/nX,'high');
    %%% pad x to mitigate edge effects
    pad = mean(x)*ones(nX,1);
    x = filtfilt(b,1,[pad;x;pad]);
    x = x(nX+1:end-nX);
end

out = struct('BFILT',[],'f',[],'dt',[],'xrec',[],'xfilt',[],'ximp',[],'a',[],'exvar',[],'B',[],'wb',[],'segment',[],'opts',opts);
xresid = x;
X = x(T0);
thresh = 1;

if opts.showprog
    fig = figure;
    im = imagesc([],tt,X);
    drawnow
end
    

Fremove = [];
skewweight = ones(size(T0,2),1);
for kk = 1:opts.ncomp
    switch opts.decomp_method
        case 'residual'
            Fremove = [];
        case 'direct'
%            xresid=x;
            Fremove=[out.f];
    end
      X = xresid(T0);
  
    if opts.resegment
       Xadj = X;
       prewin = segment.window(nX);
    else
        Xadj = diag(segment.window(nX))*X;
%        prewin = hann(nX);
      %   prewin=ones(nX,1);
         prewin=kaiser(nX,2); % Retain some tapering for stability
    end
    clear dt
    for k= 1:opts.niter + 1  %%% Apply bstd iteratively. Adding 1 because the phase response of BFILT is lagged according to the average delay in the input not output. 
       if opts.use_ideal_filter && k==opts.niter            
            [dt(k,:),Xadj,B,BFILT,wb,~,~,Bideal] = bstd(Xadj,opts.lpfilt,Fremove,prewin,skewweight,opts.bstdargs{:});
       else
           [dt(k,:),Xadj,B,BFILT,wb] = bstd(Xadj,opts.lpfilt,Fremove,prewin,skewweight,opts.bstdargs{:});
       end
       sdt = sum(dt,1);
   
       if opts.resegment 
            [T,tt] = chopper(segment.Trange, segment.wint+ round(sdt)./segment.fs,segment.fs);
            T(T<1)=1;
            T(T>n)=n;
            
            Xadj = xresid(T);
            
%             segment.wint(:,any(T<1 | T>n)) = [];
%            T(:,any(T<1 | T>n)) = [];
       end
       if opts.showprog
           set(im,'CData',Xadj*diag(skewweight))
           title(sprintf('Iter. %i',k))
           drawnow
       end
       if opts.skewness_threshold > -Inf && ~(islogical(opts.skewness_threshold)&&~opts.skewness_threshold)
           windx = round(mod(wb*size(Xadj,1),size(Xadj,1))+1);
            FXadj = fft(Xadj);
           Xfilt = real(ifft(FXadj(windx,:).*repmat(BFILT,1,size(Xadj,2))));
           skewweight = skewness(Xfilt)'>opts.skewness_threshold;
       end
       k
    end
    if opts.use_ideal_filter
        BFILT=Bideal;
    end
%     switch type
%         case 'svd'
%             [u,l] = svd(Xadj);
% 
%             [~,pckeep] = min(diff(diag(l)));
%             f = u(:,1:pckeep)*sqrt(l(1:pckeep,1:pckeep));
%         case 'mean'

            
%            f = mean(Xadj,2);
        %   sdt = sdt-median(sdt);

            sdt = sum(dt(1:end-1,:),1); %Note that BFILT accounts for the lags on the penultimate iteration, so excluding the last one here
            [T,tt] = chopper(segment.Trange, segment.wint+ round(sdt)./segment.fs,segment.fs);
            T(T<1)=1;
            T(T>n)=n;
            f =  segment.window(nX).*(xresid(T)*skewweight(:)./sum(skewweight));
%            f = Xadj*skewweight(:)./sum(skewweight);
%     end
    out(kk).f= f;

%     switch opts.decomp_method
%         case 'direct'
% %                     [T,tt] = chopper(segment.Trange, round(sum(dt,1))./segment.fs,segment.fs);
%              T(T<1)=1;
%              T(T>n)=n;
%      
%             Fremove(:,kk)= mean(xresid(T),2); %#ok<AGROW>
%     end
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
    switch opts.impulse_method
        case 'zthreshold'
            ximp = zscore(xfilt)>thresh;
        case 'kmeans'
            [km,kmc]  = kmeans(xfilt + 0./(zscore(xfilt)>1),2);
            [~,mxi] = max(kmc);
            ximp = km==mxi;
        case {'pkskew0','peak_skew0'} % Retain peaks such that remaining peaks have a 3rd cumnulant of 0;
            
            pktype = getpeak(xfilt);
            pk = find(pktype.*sign(xfilt)>0); % Keep only concave positive and convex negative peaks.
             [srt,srti] = sort(zscore(xfilt(pk)));
             pk = pk(srti);
            keepsamples = ~isnan(iterz(srt,opts.outlier_threshold)); % Suppress extreme outliers           
            m1 = cumsum(srt.*keepsamples)./cumsum(keepsamples); % cumulative mean on sorted peaks
            m2 = cumsum(srt.^2.*keepsamples)./cumsum(keepsamples); % cumulative 2nd moment
            m3 = cumsum(srt.^3.*keepsamples)./cumsum(keepsamples); % cumulative 3rd moment
            %  Third cumulant
            c3 = m3 - 3*m2.*m1 + 2*m1.^3; % Third cumulant on sorted peaks
            kept_peaks = pk(srt>0 & c3>0 ); % Keep all positive concave peaks within the set that 
         
            ximp=false(size(x));
            ximp(kept_peaks)=true;
        case 'skew0' % values such that remaining peaks have a 3rd cumnulant of 0;
            
            if opts.impulse_skewness_sd_threshold~=0
                pktype = getpeak(xfilt);
                npk = sum(pktype==1);
                skewness_threshold = opts.impulse_skewness_sd_threshold*sqrt(6/npk); % Approximate variance of skewness; use this threshold instead of 0.
            else
                skewness_threshold = 0;
            end
            
            [srt,srti] = sort(zscore(xfilt));
             keepsamples = ~isnan(iterz(srt,opts.outlier_threshold)); % Suppress extreme outliers           
            m1 = cumsum(srt.*keepsamples)./cumsum(keepsamples); % cumulative mean on sorted peaks
            m2 = cumsum(srt.^2.*keepsamples)./cumsum(keepsamples); % cumulative 2nd moment
            m3 = cumsum(srt.^3.*keepsamples)./cumsum(keepsamples); % cumulative 3rd moment
            %  Third cumulant
            c3 = m3 - 3*m2.*m1 + 2*m1.^3; % Third cumulant on sorted peaks
            keepsrt = srt>0 & c3> skewness_threshold ;
            kept_times= srti(keepsrt); % Keep all positive concave peaks within the set that 
             ximp=false(size(x));
            ximp(kept_times)=true;
%             srti(keepsrt)=[];
%             srt(keepsrt)=[];
            

    end
    
    xrec = ifft(fft(xfilt.*ximp).*fft(f));%*nX/n;
    a = xrec'*xresid./sum(xrec.^2);
    xrec = a*xrec;
   
    out(kk).BFILT = bfilt;
    %out(kk).f= mean(Xadj,2);
%     out(kk).f= Xadj*skewweight./sum(skewweight);
    
    out(kk).dt= sdt;
    out(kk).xrec = xrec;
    out(kk).xfilt = xfilt;
    out(kk).ximp= find(ximp);
    
    out(kk).a = a;
    out(kk).exvar = 1-sum(abs(xresid-xrec).^2)./sum(abs(xresid).^2);
    out(kk).B = B;
    out(kk).wb = wb;
           
    segment.wintadj = segment.wint+round(sdt)./segment.fs;
    out(kk).segment = segment;
    out(kk).segskew = skewness(Xadj);
    %%%  A simple measure of compression: 
    %%%     variance explained X total samples / number of values retained 
    %%%     (length of f + number of  impulses used in constructing xrec)
    out(kk).compression = out(kk).exvar*length(x)/(size(Xadj,1)+length(out(kk).ximp));
 %   if strcmp(decomp_method,'residual')
        xresid =xresid-xrec;
  %  end
end
    


%%%%%%%%%%
function pksign = getpeak(x)

% Get peaks and concavity

tol = eps;

dx = diff([x(end);x;x(1)]);  %Again, compute difference on circular domain

dx = dx.*(abs(dx)>tol); %Apply tolerance threshold

sdx = sign(dx) ;

dsdx = [sdx(end);sdx] - [sdx;sdx(1)];

zeroc = sign([x;x(1)])~=sign([x(end);x]);
% dsdx = diff(sdx);

pk = abs(dsdx)>0;
% pk([1 end]) = true;
pksign = sign(dsdx.*pk);

%%% When one or more adjacent values are equal, assign a peak  to
%%% the first inflection for concave or convex regions, and to both
%%% inflections for saddle points.
if any(dx==0) && any(dx~=0) 
    fpks = pksign(pk);
    pk(pk) = [fpks(end);fpks(1:end-1)] ~= fpks;
    pksign = pk.*pksign;
end
    
    

% pk = pk(2:end-1);
pksign = pksign(2:end-1);
% zeroc = zeroc(2:end-1);


