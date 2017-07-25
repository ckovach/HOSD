

function [out,Xadj,X,dt]=bsident(x,winN,lpfilt,ncomp)


if nargin < 3 || isempty(lpfilt)
    lpfilt = .1;
end
if nargin <4 || isempty(ncomp)
    ncomp = 1;
end

n = length(x);
povlp = .5;
[T,tt] = chopper([-winN/2 winN/2-1], 1:winN*(1-povlp):n,1);
T(:,[1 end]) = [];

niter = 10;

nX = size(T,1);


xresid = x;

thresh = 1;
for kk = 1:ncomp
%%
    X = xresid(T);


    Xadj = diag(hann(nX))*X;
    clear dt
    for k= 1:niter
        [dt(k,:),Xadj,B,BFILT,wb] = bstd(Xadj,lpfilt);
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
    h(n) = 0;
    h = circshift(h,-floor(nX/2));

    xfilt = ifft(fft(xresid).*fft(h));

%     pk = getpeak(xfilt);
%     pk = pk(zscore(xfilt(pk))>thresh);
% 
%     ximp = zeros(size(x));
%     ximp(pk) = 1;

    %xrec = ifft(fft(ximp).*fft(f));%*nX/n;
    xrec = ifft(fft(xfilt.*(zscore(xfilt)>thresh)).*fft(f));%*nX/n;
    a = xrec'*xresid./sum(xrec.^2);
    xrec = a*xrec;
   
    out(kk).BFILT = BFILT;
    out(kk).f= mean(Xadj,2);
    out(kk).dt= sum(dt);
    out(kk).xrec = xrec;
    out(kk).xfilt = xfilt;
    out(kk).thresh = thresh;
    
    out(kk).a = a;
    out(kk).exvar = 1-sum(abs(xresid-xrec).^2)./sum(abs(xresid).^2);
    out(kk).B = B;
    out(kk).wb = wb;
     xresid =xresid-xrec;

end
    
    
