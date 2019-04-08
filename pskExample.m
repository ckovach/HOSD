%%% This script demonstrates a simple phase-shift keying (PSK) transmission scheme.
%%% The "carrier" is a random signal within the lower portion of the spectrum and the
%%% "transmission" is a zero-phase signal within the higher portion with some specified
%%% offset relative to the carrier. The offset is less than the width of
%%% the central peak in the autocorrelation of the carrier. 
%%% 
%%% Bispectral pattern recovery is used to obtain a filter for recovering
%%% the transmission with an unknown random carrier. 

n = 1.5e5;
m = 1000; 


fs = 1e3; 


nsymb =2;
tm = ((0:m-1)'-floor(m/2))/fs;

t = ((0:n-1)'-floor(n/2))/fs;
tn = ((0:n-1)'-floor(n/2));

carrierbw = [10 20; 30 40]*2; %Lower and upper carrier bands
trbw = carrierbw([3 2]);      % Middle transmission band
% trbw = [10 20]*8    ;
% trbw2 = [36 48];

carrier = 0;
for k = 1:size(carrierbw,1)
    carrier = carrier+bpfilt(hann(m).*randn(m,1),[fs carrierbw(k,:)]);
end
% carrier = zscore(carrier); %Carrier signal
%  carrier = zscore(bpfilt(fftshift((0:m-1)'==0),[fs carrierbw]));
% carrier = zscore(bpfilt(tm==0,[fs carrierbw]));
% subn = round(1000./mean(carrierbw)/2);
% subwin = convn(tm==0,hann(subn),'same');

% trans0 = zscore(bpfilt(randn(m,1).*hann(m).^4,[fs trbw]));
% trans1 = zscore(bpfilt(randn(m,1).*hann(m).^4,[fs trbw]));
%  trans0= zscore(bpfilt(diff(fftshift((0:m)'==0)),[fs trbw]));
carrier_center_frequency = mean(carrierbw);
peakw = 1./carrier_center_frequency/8;
trans = [];
% tr = zscore(bpfilt(randn(m,1).*hann(m),[fs trbw]));
tr = bpfilt(randn(m,1).*hann(m),[fs trbw]);
% tr = zscore(bpfilt((tm==0).*hann(m),[fs trbw]));

%%% Create the phase-shifte transmission symbols
trph = linspace(-pi/2,pi/2,nsymb+1); % Use the half-range
trans=[];
transfft = fft(tr);
for k = 1:nsymb
%    trans(:,k) = zscore(bpfilt(tm==round(1000.*peakw*(1-2*(k-1)/(nsymb-1)))/1000,[fs trbw]));
%     trans(:,k) = zscore(convn(tr,tm==round(1000.*peakw*(1-2*(k-1)/(nsymb-1)))/1000,'same'));
%     trans(:,k) = zscore(convn(tr,tm==0,'same'));
     trans(:,k)= real(ifft(transfft.*exp(sign(fftshift(tm))*1i*trph(k+1))));
end
% 
% cph = linspace(-pi/2,pi/2,nsymb);
 carr=[];
% carrfft = fft(carrier);

for k = 1:nsymb
%     carr(:,k) =real(ifft(carrfft.*exp(sign(fftshift(tm))*1i*cph(k))));
    carr(:,k) = carrier;
end
    % trans0 = zscore(bpfilt(tm==-round(1000.*peakw)/1000,[fs trbw]));
%  trans1 = -trans0;
% trans1 = zscore(bpfilt(subwin.*randn(m,1),[fs trbw]))/sqrt(m./subn);
% trans0 = zscore(bpfilt(subwin.*randn(m,1),[fs trbw]))/sqrt(m./subn);
% 
% B0 = carrier + trans0;
% B1 = carrier + trans1;

% B = repmat(carrier,1,nsymb)+trans;
B = zscore(carr+trans);

TRt = zeros(n,1);
rate = 2;
dtrt = ones(round(rate*n/m),1)*m./rate;
dtrt = dtrt+round(randn(length(dtrt),1)*m/8);
trt = round(cumsum(dtrt));
trt(trt>n)=[];
TRt(trt)=1;
TRt(end+1:n)=0;
trsym = TRt.*ceil(rand(n,1)*nsymb);
N = zeros(n,nsymb);
N(find(TRt) + n*(trsym(find(TRt))-1))=1;
% N = rand(n,nsymb)<2/(1*m)/nsymb;
Ssep = [];
for k = 1:nsymb
    Ssep(:,k) = convn(N(:,k),B(:,k),'same');
end
 % N0 = rand(n,1)<1/m;
% N1 = rand(n,1)<1/m;

% S0 = convn(N0,B0,'same');
% S1 = convn(N1,B1,'same');
% S = S0+S1;
% S = zscore(sum(Ssep,2));
S =sum(Ssep,2);

%%% Multipath fading
% S = S+.5*circshift(S,round(m*.25));


arg = @(x)atan2(imag(x),real(x))

noiseamp = -6; %Signal to noise ratio in dB
GN =10^(-noiseamp/20)*zscore(lpfilt(randn(n,1),[fs 1.25*max([trbw(:);carrierbw(:)])*2]));

X = S+GN;

ncomp =1;
clear hos
hos(ncomp) = hosobject;
hos.initialize(m,fs,max([trbw(:);carrierbw(:)])*1.25);
[Xsh,Xwin] = hos.get_block(X);


%%
xrec = hos.xrec(X);
xfilt = hos.xfilt(X);
ximp = full(hos.ximp(X));
xthresh = full(hos.xthresh(X));
%%% Time spread of the feature
% ffun = ifftshift(real(ifft(hos(1).filterfft.*hos(1).wavefft)));
ffun = hos(1).filterfun;
ftw = sqrt(tm'.^2*[ffun].^2./sum([ffun].^2) -(tm'*[ffun].^2./sum([ffun].^2)).^2 );

%%% Smooth according to the timespread
ff = hann(round(mean(ftw).*fs*2));
ff=ff./sqrt(sum(ff));
xthrf = convn(xthresh, ff,'same');
[pk,pks] = getpeak2(xthrf(:,1));
pkt = find(pks>0)/fs;
[T,tt] = chopper([-1 1]*m/fs/2,pkt,fs);
T(T<1)=1;
T(T>size(xrec,1))=size(xrec,1);
subT = T(abs(tt)<.1,:);

count = arrayfun(@(k)sum(N(subT + (k-1)*size(N,1))),1:nsymb,'uniformoutput',false);
count = cat(1,count{:});
[mx,type] = max(count);
type(mx==0)=0;
[srt,srti] = sort(type);

% xfilt0 = hos.xfilt(S0+GN);
% xfilt1 = hos.xfilt(S1+GN);
% ximp0 = hos.ximp(S0+GN);
% ximp1 = hos.ximp(S1+GN);

% xff = bpfilt(xfilt(:,1),[fs trbw]);
% hxff = hilbert(xff);
% h0 = hist(arg(hxff(ximp0(:,1))),-pi:.1:pi);
% h1 = hist(arg(hxff(ximp1(:,1))),-pi:.1:pi);

% hxff = hilbert(bpfilt(xfilt(:,1),[fs carrierbw]));
hxff = hilbert(bpfilt(xfilt(:,1),[fs trbw]));

h=[];
ph={};
th=-pi:.1:pi;
for k = 1:nsymb
    ph{k} = hxff(round(pkt(type==k)*fs));
    h(:,k) = hist(arg(hxff(round(pkt(type==k)*fs))),th);
% h2 = hist(arg(hxff(round(pkt(type==2)*fs))),th);
end

wb = fftshift(hos(1).freqindx.Bfreqs{1});
figure
subplot(2,2,1)
% plot(tm,B)
plot(tm,trans + carr - 10*std(trans(:)+carr(:))*ones(size(tm))*(-ceil(nsymb/2)+1:floor(nsymb/2)))
hold on
plot(tm,carr - 10*std(trans(:)+carr(:))*ones(size(tm))*(-ceil(nsymb/2)+1:floor(nsymb/2)),'k')
title('Transmission Signals')
legend([arrayfun(@(k)sprintf('Symbol %i',k),1:nsymb,'uniformoutput',false),{'Carrier'}])

subplot(2,2,2)
imagesc(wb,wb,fftshift(abs(hos(1).bicoh)))
axis xy image
title('Bicoherence')

subplot(2,2,3)
% imagesc(tm,[],xfilt(T(:,srti))')
 imagesc(tm,[],real(hxff(T(:,srti)))')
title(sprintf('Filtered vs Symbol, SNR: %0.0f dB',noiseamp))
set(gca,'ytick',1:4:length(srt),'yticklabel',srt(1:4:end))
hold on
plot([0 0],[0 size(T,2)],'k')
ylabel 'Symbol'
xlabel 'time (s)'
xlim([-1 1]*.1)

subplot(2,4,8)
% stem(th,h)
% xlim([-pi pi]/2)
% xlabel 'Phase (rad)'
% title('Distribution of Phase in the Carrier band')
% % legend({'Symbol 1','Symbol 2'})
% legend(arrayfun(@(k)sprintf('Symbol %i',k),1:nsymb,'uniformoutput',false))
for k = 1:nsymb
    hold on
    plot(hxff(round(pkt(type==k)*fs)),'.')
end
% xlim([-pi pi]/2)
% xlabel 'Phase (rad)'
title('Distribution of Phase in the Transmission band')
% legend({'Symbol 1','Symbol 2'})
legend(arrayfun(@(k)sprintf('Symbol %i',k),1:nsymb,'uniformoutput',false))


[u,l,v] = svd(detrend(xfilt(T)',0)');

subplot(2,4,7)
title('PCA')
for k = 1:nsymb
    hold on
    plot(v(type==k,1),v(type==k,2),'.')
end
xlabel PC1
ylabel PC2
% legend(arrayfun(@(k)sprintf('Symbol %i',k),1:nsymb,'uniformoutput',false))


p = h./repmat(sum(h),size(h,1),1);
dp = repmat(p,1,size(p,2)).*log(repmat(p+eps,1,size(p,2))./kron(p+eps,ones(1,size(p,2))));
kldist  = reshape(sum(dp),[1 1]*nsymb);
