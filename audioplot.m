
ncascade = 4;

h = dsp.AudioRecorder;
h.SampleRate = 8000;
h.SamplesPerFrame = 1200;
h.NumChannels=1;
Fs = h.SampleRate;

for k = 1:ncascade
hos = hosobject(3,600,h.SampleRate,h.SampleRate/2);
%hos.window= @(x)kaiser(x,3);

hos.glowpass = h.SampleRate/2*.99;
hos.hos_learning_rate = .001;
hos.burnin = 1./hos.hos_learning_rate;
hos.filter_adaptation_rate = .001;
% Xn = ecgfn(T);
hos.reset
hos.update_frequency_indexing
hoss(k)=hos;
end
X = step(h);

tt = (0:length(X)-1)/Fs;

figure, 
subplot(2,2,1), 
w = fftshift(hos.freqs{1});
kpw = false(size(w));
kpw(hos.keepfreqs{1}) = true;
w = w(fftshift(kpw));
im(1) = imagesc(w,w,fftshift(abs(hos.Bfull(:,:)))); 
%caxis([0 1])
colorbar
axis image
subplot(2,2,2), 
 imwf = plot(ifftshift(hos.sampt)/Fs,zscore([hoss.waveform]) + 2*ones(hos.buffersize,1)*[1:ncascade]); 
% im(2:3) = plot(ifftshift(hos.sampt)/Fs,[hos.waveform,hos.shiftbuffer]);
% set(im(2),'color',[1 1 1]*.75)
axis tight
% ylim([-1 1]*10)
grid on
subplot(2,1,2), 
im(4:5) = plot(tt,X(:,[1 1 ])); 
set(im(4),'color',[1 1 1]*.75)
% set(im(5),'color',[1 1 1]*0,'linewidth',2)
set(im(5),'linewidth',2,'color','r')
grid on
axis tight
ylim([-1 1]*.5)

nfr = 35;
%%
%   Xs = zeros(h.SamplesPerFrame,nfr);
Xrs = permute(Xs(:,:,ones(1,ncascade)),[1 3 2]);
Xfs = Xrs;
Xresids = Xs;
go = true;
clear Xr
for kkk = 1:nfr*1e6
    k = mod(kkk-1,nfr)+1
%       Xn = step(h);
   Xn = Xs(:,k);
     Xresid = Xn;
    
    for kk = 1:ncascade
        hoss(kk).get_input(Xresid);
        Xr(:,kk) = hoss(kk).reconstruct(Xn);
        Xf(:,kk) = hoss(kk).apply_filter(Xn);
        Xresid = Xresid-Xr(:,kk);
    end
    
    set(im(1),'CData',fftshift(abs(hoss(1).bicoh(:,:)))), 
   % set(im(2), 'Ydata',fftshift(hoss(1).shiftbuffer)*0), 
%    set(im(3), 'Ydata',fftshift([hoss.waveform]));
    for kk = 1:ncascade
        set(imwf(kk), 'Ydata',fftshift(zscore([hoss(kk).waveform]),1) + 5*kk);
    end
    set(im(4), 'Ydata',Xn);
%     set(im(5), 'Ydata',X(:,k));
 
     set(im(5), 'Ydata',sum(Xr,2));
    %k%[hos.EDF hos.current_threshold]
%     [k, hos.EDF],
%     pause(.1),
    drawnow,
    Xs(:,k) = Xn;
    Xrs(:,:,k) = Xr;
    Xresids(:,k) = Xresid;
    Xfs(:,:,k) = Xf;
    %Xfs(:,k)= hos.apply_filter(Xn);
end

xrs = reshape(permute(Xrs,[1 3 2]),size(Xrs,1)*size(Xrs,3),size(Xrs,2));

