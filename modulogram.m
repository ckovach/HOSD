function [MGNRM,bfreq,mfreq,out] = modulogram(x,Nwin,Fs,Novlp,lowpass)


if nargin < 3 || isempty(Fs)
    Fs = 1;
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
window = window/sum(window);
x = x-nanmean(x);

T = chopper([0 Nwin-1],0:Novlp:length(x)-Nwin,1,length(x));
X = x(T).*window;

disc = any(isnan(X));
if any(disc)
    fprintf('\n%i segments with nan values discarded (%0.2f%%)',sum(disc),mean(disc)*100);
end
X = X(:,~disc);

X = [X;zeros(size(X))];

w = ifftshift((0:2*Nwin-1)-floor(2*Nwin/2));

[MF,BF] = ndgrid(w(w>=0 & w <=lowpass(1)),w(w>=0 & w <= lowpass(2)));

W1 = BF(1:2:end,:);
W2 = -BF(1:2:end,:) + MF(1:2:end,:);
W3 = -BF(1:2:end,:) - MF(1:2:end,:);

FX = fft(X);

Ws = {W1,W1,W2,W3};

PSD = mean(abs(FX).^2,2);

winFT = fft([window;zeros(size(window))]);

% DR1 = W1==0 | W2==0 | W3==0;
DR2 = W1+W2==0 | W1+W3==0 | W2+W3==0 | W1==0 ;

MGX = 1;
PSDX = 1;
winHOS= 1;
xx = 0;
for k = 1:length(Ws)
   
   I = mod(Ws{k},2*Nwin)+1;

   MGX =  MGX.*FX(I(:),:);
    
   PSDX = PSDX.*PSD(I);
   
   winHOS = winHOS.*winFT(I);
   xx = xx+w(I);
end
 
% 
% DGcorr = ifft2(fft2(abs(winHOS)).*fft2(DR2.*sqrt(PSDX)));

 MGX(DR2,:) = MGX(DR2,:) - sqrt(PSDX(DR2));

MG = mean(MGX,2);% - DGcorr(:);
NRM = mean(abs(MGX),2);
% NRM = sqrt(PSDX(:));
MGNRM = MG./NRM;

if nargout >3
   [b,dev,pval,iXX,sigma,res,Yfit,df] = complexglm(MGX(:,1:2:end)',[],'diagonly',false);
   out.beta = reshape(b,size(I));
   out.pval = reshape(pval,size(I));
   out.iXX  = iXX;
   out.se = reshape(iXX*sigma ,size(I));
   out.MG = reshape(MG,size(I));
   out.NRM = reshape(NRM,size(I));
end

% MG = reshape(MG,size(I));
MGNRM = reshape(MGNRM,size(I));

% bfreq =BF(1,1:2:end)/(2*Nwin)*Fs;
bfreq =BF(1,:)/(2*Nwin)*Fs;
mfreq =MF(1:2:end,1)/(2*Nwin)*Fs;







