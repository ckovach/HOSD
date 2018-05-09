function [dt,Xadj,Bout,Bfilt,w,NRM,BIAS] = bstd(X,lowpass,Fremove,prewin,win_weight,postwin,highpass)

% Bispectral time delay estimation.
%
% [dt,Xadj,B,BFILT] = bstd(X,lowpass)
%
%  Inputs:
%       X : data matrix in time x trials  form
%       lowpass : lowpass cutoff as a proportion of sampling rate. (NB 
%                 proportion of sampling rate, not nyquist)
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
remove_supernyq = true; % Filter bispectrum in super-nyquist region.
persistent A wlast

n = size(X,1);
m = size(X,2);
%fwin = rectwin(n);
%fwin = hann(n);
if nargin < 5 || isempty(prewin)
    prewin = kaiser(n,3);  %%% Default windowing in the time domain prior to calculating the bispectum
end
sdpostwin = .25;
if nargin < 6 || isempty(postwin)
    postwin = @(x)exp(-(x./sdpostwin).^2); %%% The bispectrum is smoothed according to this time-domain window (applied to the 3rd moment)
elseif isnumeric(postwin)
    postwin = @(x)exp(-(x./postwin).^2);
end

if nargin < 2
    lowpass = .25;
end
if nargin < 6 || isempty(win_weight) || isscalar(win_weight)
   win_weight = ones(size(X,2),1); 
else
    win_weight = win_weight(:);
end
wfull = ifftshift((0:n-1) - floor(n/2))/n;
w = wfull(abs(wfull)<=lowpass);
nb = length(w);
%twin = fftshift(hann(nb));
if nargin <7
    highpass=w(4);
end

t = (fftshift((0:length(w)-1)-ceil(length(w)/2)))./length(w);
[W1,W2] =ndgrid(w,w);

PDConj = (W1<0 | W2<0) & W1+W2>0 | W1<=0&W2<=0 ; %%% Conjugate symmetric part

PDIndx = (W1>=0 & W2>=0 | W1<=0 & W2>=0 & W1+W2<=0) & ~PDConj; %%% The cross bispectrum is symmetric with 
                                  %%%respect to only the exchange of w2 and -w2-w1, 
                                  %%%so need to include the entire positive
                                  %%%quadrant and one region in W1<0 & W2>0
                                  %%% quadrant.
%%% Map into principal domain
W1pd=W1;W2pd=W2;
W1pd(PDConj) = -W1(PDConj);
W2pd(PDConj) = -W2(PDConj);
wpi = W1pd >= 0 & W2pd<=0 & W1pd + W2pd <=0; 
W2pd(wpi) = -W1pd(wpi)-W2pd(wpi);

%%% Indices into the bispectral plane
I1pd = round(mod(W1pd*n,nb)+1);
I2pd = round(mod(W2pd*n,nb)+1);

% Map of all regions into the principal domain
pdmap = sub2ind([1 1]*nb,I1pd,I2pd);



[T1,T2] =ndgrid(t,t);

%kpi = abs(W1)<lowpass & abs(W2)<lowpass;

W3 = mod(-W1-W2+.5,1)-.5;
supernyq = abs(W1+W2)>.5; %region where frequency sum exceeds Nyquist.
I1 = round(mod(W1(PDIndx)*n,n)+1);
I2 = round(mod(W2(PDIndx)*n,n)+1);
I3 = round(mod(W3(PDIndx)*n,n)+1);

windx = round(mod(w*n,n)+1);


FX = fft(X);
FXwin = fft(X.*repmat(prewin(:),1,m));
BB = FXwin(I1,:).*FXwin(I2,:).*FXwin(I3,:); 

%%% Average Bispectrum
NRM = zeros(nb);
BIAS = NRM;
switch normalization
    case 'awplv' %%% This normalizes such that the result is an amplitude-weighted mean phase difference. See Kovach 2017 (IEEE Trans. Sig. Proc. v65 n17, p4468)
%        nrm= sum(abs(BB),2);
%        BIAS(PDIndx) = sqrt(sum(abs(BB).^2,2)./(nrm.^2+eps));
        nrm = abs(BB)*win_weight(:);
        BIAS(PDIndx) = sqrt((abs(BB).^2*win_weight)./(nrm.^2+eps));
        BIAS(:) = BIAS(pdmap);
    case 'bicoh' %%% Normalize according to standard bicoherence
        BIAS=0;
%         nrm = sqrt(sum(abs(FXwin(I1(:),:)).^2,2).*sum(abs(FXwin(I2(:),:).*FXwin(I3(:),:)).^2,2));
         nrm = sqrt((abs(FXwin(I1(:),:)).^2*win_weight).*(abs(FXwin(I2(:),:).*FXwin(I3(:),:)).^2*win_weight));
%         NRM(:) = sqrt(sum(abs(FX(I3(:),:)).^2,2).*sum(abs(FX(I2(:),:).*FX(I1(:),:)).^2,2));
    case 'none'  %%% No normalization.
        BIAS=0;
        nrm = size(BB,3);
    otherwise
        error('unrecognized normalization')
end
NRM(PDIndx) = nrm;
NRM = NRM(pdmap);

% sBBnrm = sum(BB,2)./nrm;
sBBnrm = (BB*win_weight)./nrm;
if nargin >2 && ~isempty(Fremove)
    
    FFremove = fft(Fremove.*repmat(prewin(:),1,size(Fremove,2)));
%    FFremove = fft(Fremove);
    Bremove =  FFremove(I1,:).*FFremove(I2,:).*FFremove(I3,:)./repmat(nrm,1,size(FFremove,2)); 
    Bremove = Bremove*(Bremove'*Bremove)^-.5;
    sBBnrm = sBBnrm- Bremove*(Bremove'*sBBnrm);
end
B = zeros(nb);
B(PDIndx) = sBBnrm;
B = B(pdmap);
B(PDConj) = conj(B(PDConj));

B = B-BIAS.*B./(abs(B)+eps);
B(isnan(B))=0;
Bout=B;


if snr_weighting %%% This option attempts to estimate the signal-to-noise ratio and adjust weighting accordingly
                 %%% More computationaly intensive, and doesn't seem to
                 %%% produce appreciably improved results.
    
    if isempty(A)|| ~isequal(wlast,wfull)
      A = sparse(size(X,1),sum(PDIndx(:)),3*sum(PDIndx(:)));
        Aindx = mod(wfull([I1(:),I2(:),I3(:)])*size(X,1),size(X,1)) + 1 +repmat(((1:sum(PDIndx(:)))'-1)*size(X,1),1,3);
        A(Aindx(:))=1;
        wlast = wfull;
    end
    TMW=zeros(size(B));
    TMW(:) = postwin(T1).*postwin(T2).*postwin(T1-T2);
    B = fft2(ifft2(B).*TMW);
    B(abs(W1)<highpass|abs(W3)<highpass|abs(W2)<highpass| abs(W3)>=lowpass | (remove_supernyq & supernyq))=0;
  
    fls = zeros(size(X,1),1);
    geti = abs(wfull)<=max(abs(w));
    fls(geti)=  log(abs(B(PDIndx)).^2+eps)'/A(geti,:);
   % wfls = wfull(abs(wfull)<=max(abs(w))*2);
    fls(abs(wfull)>max(abs(w)))=log(eps);
    snr = exp(fls)./(1-exp(fls));
 %   SNR = snr(I1).*snr(I2);%.*snr(I3);
    SNR = zeros(size(B));
    SNR(PDIndx) = snr(I1).*snr(I2).*snr(I3);
    SNR = SNR(pdmap);
    
 %   nrm = @(x)x./(abs(x)+eps).*abs(B).^2./(1-abs(B).^2);
    nrmfun = @(x)x./(abs(x)+eps).*SNR;
else    
     B(:) = B(:)./(NRM(:)+eps);
     B(abs(W1)<highpass|abs(W3)<highpass|abs(W2)<highpass | abs(W3)>=lowpass | (remove_supernyq & supernyq))=0;
  
     TMW=zeros(size(B));
     TMW(:) = postwin(T1).*postwin(T2).*postwin(T1-T2);
     B = fft2(ifft2(B).*TMW);
    nrmfun = @(x)x;
%    nrm = @(x)x./(abs(x)+eps).*abs(B).^2;
    
end

B23 = zeros(size(B));
% B23(PDIndx) = sum(FXwin(I2(:),:).*FXwin(I3(:),:),2);
B23(PDIndx) = (FXwin(I2(:),:).*FXwin(I3(:),:)*win_weight);
B23 = B23(pdmap);
B23(PDConj) = conj(B23(PDConj));

%%% The optimal filter
Bfilt = sum(nrmfun(B23.*conj(B)),2);
Bfilt(isnan(Bfilt))=0;

%%% Delays are determined from maxima in the data filtered by Bfilt.
FPH = ifft(repmat(Bfilt(:),1,m).*FXwin(windx,:));
[~,mxi] =max(real(FPH));
dt = w(mxi)./max(w)*n/2;
dt = dt-mean(dt);

%%% Circular shift according to the delay from the previous step
FXadj = FX.*exp(1i*wfull'*dt*2*pi);
Xadj = real(ifft(FXadj));


