function [dt,Xadj,Bout,Bfilt,w,NRM,BIAS] = bstd(X,lowpass,tmwin,highpass)

% Bispectral time delay estimation.
%
% [dt,Xadj,B,BFILT] = bstd(X,lowpass)
%
%  Inputs:
%       X : data matrix in time x trials  form
%       lowpass :   lowpass cutoff as a proportion of sampling freq.
%   
%  Outputs:
%       dt : Estimated relative time delays for each trial
%       Xadj : Data shifted according to dt
%       B : Signal bispectum
%       BFILT: Optimal filter function for feature extraction
%
% C. Kovach 2017

normalization  = 'awplv';
snr_weighting = false;

persistent A wlast

n = size(X,1);
m = size(X,2);
%fwin = rectwin(n);
%fwin = hann(n);
fwin = kaiser(n,3);
sdtmwin = .25;
if nargin < 3 || isempty(tmwin)
    tmwin = @(x)exp(-(x./sdtmwin).^2);
elseif isnumeric(tmwin)
    tmwin = @(x)exp(-(x./tmwin).^2);
end

if nargin < 2
    lowpass = .25;
end

wfull = ifftshift((0:n-1) - floor(n/2))/n;
w = wfull(abs(wfull)<=lowpass);
nb = length(w);

if nargin <4
    highpass=w(4);
end

%twin = fftshift(hann(nb));

t = (fftshift((0:length(w)-1)-ceil(length(w)/2)))./length(w);
[W1,W2] =ndgrid(w,w);
[T1,T2] =ndgrid(t,t);

%kpi = abs(W1)<lowpass & abs(W2)<lowpass;

W3 = mod(-W1-W2+.5,1)-.5;

I1 = round(mod(W1*n,n)+1);
I2 = round(mod(W2*n,n)+1);
I3 = round(mod(W3*n,n)+1);

FX = fft(X);
FXwin = fft(X.*repmat(fwin(:),1,m));
BB = FXwin(I1(:),:).*FXwin(I2(:),:).*FXwin(I3(:),:); 

%%% Average Bispectrum
NRM = zeros(nb);
BIAS = NRM;
switch normalization
    case 'awplv'
        NRM(:) = sum(abs(BB),2);
        BIAS(:) = sqrt(sum(abs(BB).^2,2)./NRM(:).^2);
    case 'bicoh'
        BIAS=0;
         NRM(:) = sqrt(sum(abs(FXwin(I1(:),:)).^2,2).*sum(abs(FXwin(I2(:),:).*FXwin(I3(:),:)).^2,2));
%         NRM(:) = sqrt(sum(abs(FX(I3(:),:)).^2,2).*sum(abs(FX(I2(:),:).*FX(I1(:),:)).^2,2));
    case 'none'
        BIAS=0;
        NRM = size(BB,3);
    otherwise
        error('unrecognized normalization')
end

B = zeros(nb);
B(:) = sum(BB,2)./NRM(:);

B = B-BIAS.*B./abs(B);
Bout=B;


if snr_weighting
    
    
    if isempty(A)|| ~isequal(wlast,wfull)
      A = sparse(size(X,1),numel(B),3*numel(B));
        Aindx = mod(wfull([I1(:),I2(:),I3(:)])*size(X,1),size(X,1)) + 1 +repmat(((1:numel(B))'-1)*size(X,1),1,3);
        A(Aindx(:))=1;
        wlast = wfull;
    end
    TMW=zeros(size(B));
    TMW(:) = tmwin(T1).*tmwin(T2).*tmwin(T1-T2);
    B = fft2(ifft2(B).*TMW);
    B(abs(W1)<highpass|abs(W3)<highpass|abs(W2)<highpass)=0;
  
    fls = zeros(size(X,1),1);
    geti = abs(wfull)<=max(abs(w));
    fls(geti)=  log(abs(B(:)).^2+eps)'/A(geti,:);
   % wfls = wfull(abs(wfull)<=max(abs(w))*2);
    fls(abs(wfull)>max(abs(w)))=log(eps);
    snr = exp(fls)./(1-exp(fls));
 %   SNR = snr(I1).*snr(I2);%.*snr(I3);
    SNR = snr(I1).*snr(I2).*snr(I3);
    
 %   nrm = @(x)x./(abs(x)+eps).*abs(B).^2./(1-abs(B).^2);
    nrm = @(x)x./(abs(x)+eps).*SNR;
else    
     B(:) = B(:)./(NRM(:));

     TMW=zeros(size(B));
     TMW(:) = tmwin(T1).*tmwin(T2).*tmwin(T1-T2);
     B = fft2(ifft2(B).*TMW);
     B(abs(W1)<highpass|abs(W3)<highpass|abs(W2)<highpass | abs(W3)>=lowpass)=0;
    nrm = @(x)x;
%    nrm = @(x)x./(abs(x)+eps).*abs(B).^2;
    
end

B23 = reshape(sum(FXwin(I2(:),:).*FXwin(I3(:),:),2),size(B));

Bfilt = sum(nrm(B23.*conj(B)),2);
FPH = ifft(repmat(Bfilt(:),1,m).*FXwin(I1(:,1),:));
[~,mxi] =max(real(FPH));
dt = w(mxi)./max(w)*n/2;
dt = dt-mean(dt);

FXadj = FX.*exp(1i*wfull'*dt*2*pi);
Xadj = real(ifft(FXadj));


