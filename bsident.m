

function [out,Xadj,X,dt]=bsident(x,segment,lpfilt,ncomp,varargin)

niter = 10;

%%% Method to identify impulse
%impulse_method = 'zthreshold';
impulse_method = 'kmeans';


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


xresid = x;

thresh = 1;
for kk = 1:ncomp
%%
    X = xresid(T);


    Xadj = diag(segment.window(nX))*X;
    clear dt
    for k= 1:niter
        [dt(k,:),Xadj,B,BFILT,wb] = bstd(Xadj,lpfilt,varargin{:});
        k
    end

    % FX = fft(X);
    % wfull = ifftshift((0:nX-1) - floor(nX/2))/nX;
    % FXshift = FX.*exp(1i*wfull'*sum(dt)*2*pi);
    % Xshift = real(ifft(FXshift));
    % f = mean(Xshift,2);
    f = mean(Xadj,2);

    f(n) = 0;
    f = circshift(f,-floor(nX/2));

    H = fftshift(BFILT);
    %H = H./sqrt(sum(abs(H).^2));
    H(nX) = 0;
    H = circshift(H,-floor(length(BFILT)/2));
    h = real(ifft(H));
    bfilt = h;
    h(n) = 0;
    h = circshift(h,-floor(nX/2));

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
    segment.wintadj = round(segment.wint+sum(dt));
    out(kk).segment = segment;

     xresid =xresid-xrec;

end
    
    
