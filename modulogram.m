function [MGNRM,bfreq,mfreq,out,MGX,X] = modulogram(x,Nwin,Fs,Novlp,lowpass,NFFT,highpass, NORMX,compute_bias)


if nargin < 3 || isempty(Fs)
    Fs = 1;
end

if size(x,2)>1
    Nwin = size(x,1);
    NFFT = size(x,1);
end

if isscalar(Nwin)
    window = rectwin(Nwin);
else
    window = Nwin;
    Nwin = length(window);
end

if nargin < 5
    lowpass = Fs/2*[1 1];
elseif isscalar(lowpass)
    lowpass = [1 1]*lowpass;
end
if nargin < 4 || isempty(Novlp)
    Novlp = ceil(Nwin*.5);
elseif Novlp<1
    Novlp = ceil(Nwin*Novlp);
end
if nargin < 6
    NFFT = Nwin;
end
if nargin < 7
    highpass = [0 0];
elseif isscalar(highpass)
    highpass = [1 1]*highpass;
end

if nargin < 8 || isempty(compute_bias)
    compute_bias = true;
end

w = ifftshift((0:NFFT-1)-floor(NFFT/2));

window = window/sum(window);
x = x-nanmean(x(:));

if size(x,2) >1
    X = x.*window;
else
    T = chopper([0 Nwin-1],0:Novlp:length(x)-Nwin,1,length(x));
    X = x(T).*window;
end

X(end+1:NFFT,:) = 0;


disc = any(isnan(X));
if any(disc)
    fprintf('\n%i segments with nan values discarded (%0.2f%%)',sum(disc),mean(disc)*100);
end
X = X(:,~disc);
if size(x,2)==1
    T = T(:,~disc);
end

FX = fft(X);


if nargin < 8 || isempty(NORMX)    
%     NORMFX = FX;
    NORMFX=[];
elseif size(NORMX,2)==1
     T2 = chopper([0 Nwin-1],0:Novlp:length(NORMX)-Nwin,1,length(NORMX));
    NORMX = NORMX(T2).*window;
    NORMX(end+1:NFFT,:) = 0;   
    NORMFX = fft(NORMX(T2));
else
    NORMFX = fft(NORMX);
end
    
% X = [X;zeros(size(X))];


nlowpass = lowpass*size(X,1)/Fs;
nhighpass = highpass*size(X,1)/Fs;

[MF,BF] = ndgrid(w(w>=nhighpass(1) & w <=nlowpass(1)),w(w>=nhighpass(2) & w <= nlowpass(2)));

W1 = BF;
W2 = -BF + MF;
W3 = -BF - MF;


Ws = {W1,W1,W2,W3};

PSD = mean(abs(FX).^2,2);

padwin = window;
padwin(end+1:NFFT) = 0;
winFT = fft([padwin;zeros(size(padwin))]);

% DR1 = W1==0 | W2==0 | W3==0;
DR2 = W1+W2==0 | W1+W3==0 | W2+W3==0 | W1==0 ;

MGX = 1;
PSDX = 1;
winHOS= 1;
NORMGX=1;
xx = 0;
for k = 1:length(Ws)
   
   I = mod(Ws{k},NFFT)+1;

   MGX =  MGX.*FX(I(:),:);
    if ~isempty(NORMFX)    
        NORMGX =  NORMGX.*NORMFX(I(:),:);
    end
   PSDX = PSDX.*PSD(I);
   
%    winHOS = winHOS.*winFT(I);
   xx = xx+w(I);
   II(:,k) = I(:);
   
end
if isempty(NORMFX)
    NORMGX = MGX;
end
% 
% DGcorr = ifft2(fft2(abs(winHOS)).*fft2(DR2.*sqrt(PSDX)));

 MGX(DR2,:) = MGX(DR2,:) - sqrt(PSDX(DR2));

MG = mean(MGX,2);% - DGcorr(:);
NRM = mean(abs(NORMGX),2);
% NRM = sqrt(PSDX(:));
BIAS = compute_bias*sum(abs(NORMGX).^2,2)./sum(abs(NORMGX),2).^2;
MGNRM = MG./NRM;
MGNRM = MGNRM./abs(MGNRM).*(abs(MGNRM)-sqrt(BIAS));

if nargout >3
   [b,dev,pval,iXX,sigma,res,Yfit,df,sigmarel] = complexglm(MGX',[],'diagonly',false);
   out.beta = reshape(b,size(I));
   out.pval = reshape(pval,size(I));
   out.iXX  = iXX;
   out.se = reshape(iXX*sigma ,size(I));
   out.MG = reshape(MG,size(I));
   out.NRM = reshape(NRM,size(I));
   out.serel = reshape(iXX*sigmarel,size(I)); %To find variance of real and imaginary parts, Re{out.se+out.srel}/2 and Re{out.se - out.srel}/2, respectively.
end

% MG = reshape(MG,size(I));
MGNRM = reshape(MGNRM,size(I));

% bfreq =BF(1,1:2:end)/(2*Nwin)*Fs;
bfreq =BF(1,:)/(NFFT)*Fs;
mfreq =MF(:,1)/(NFFT)*Fs;
% 
% if nargout >5
%    %%% Efficient partial delay filters for the diagonal slice (not working yet) 
%    fPDfilt = zeros(size(X));
%    ZFX = FX./sqrt(PSD+eps + mean(PSD));
%    for k = 1:size(X,2)
%        XXF = ZFX(:,k).*conj(ZFX);
%        fPDfilt(:,k) = mean(XXF.*conj(ZFX).*fft(ifft(XXF).^2),2);
%    end
%    t = ifftshift((0:NFFT-1)' - floor(NFFT/2));
% %    FXsh = FX;
%    PDFsh = fPDfilt;
%    for k = 1:25
%       mfpdf = mean(PDFsh,2);
%       Xdelay = ifft( mfpdf.*FX);
%       [~,mxi] = max(abs(Xdelay));
%       delt = t(mxi)';
% %        delt = delt - round(angle(mean(exp(1i*2*pi*delt/NFFT)))/(2*pi)*NFFT);
%       D = t==delt;
% %       D = sign(D.*Xdelay);
%       FD = fft(D);
% %       FXsh = FXsh.*FD;
%       PDFsh = fPDfilt.*conj(FD).^2.*sign(sum(D.*Xdelay));
%    end
%    FD = FD*(-1)^(sum(D(:))<0);
%    Xsh = ifft(FX.*conj(FD));
%    mpdf = ifft(mfpdf);
%    [~,mxi] = max(abs(ifft(abs(mfpdf).*mean(FX.*conj(FD),2))));
%    mpdf = fftshift(circshift(mpdf,t(mxi)),1);
%    Xsh = fftshift(circshift(Xsh,-t(mxi)),1);
% end
% 
% 



