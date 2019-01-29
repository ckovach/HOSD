

% This script demonstrates HOSD3 on records retrieved from the PhysionetDB
% database. It expects the WFDB Toolbox for Matlab and Octave, available at
%    https://www.physionet.org/physiotools/matlab/wfdb-app-matlab/
%
% Examples are based on those in Kovach and Howard 2019


d = input('Which demo to run?\n1) Normal sinus rhythm in Gaussian noise.\n2) Normal sinus rhythm in non-Gaussian noise.\n3) Abnormal rhythm.\n[1]:');

demo_list = {'nsr_ecg_noise'  % Denoising of Normal ECG
            'nsr_ecg_chi2_noise'  % Denoising of Normal ECG with non-Gaussian (chi-square) noise
            'arrhythmia'};  % Component separation of abnormal ECG

if isempty(d), d=1;end
example = demo_list{d};

window_duration = 1.5; %Window duration in seconds.
    
clim =[0 .5];
lowpass= 64;
switch example
    case 'nsr_ecg_noise'  % Denoising of Normal ECG
        
        db = 'nsrdb'; %% Database to query
        recordID =   '16483'; %% Record ID to get
        InbandNoise = 4; % In-band noise in dB
        OutbandNoise = 7; %Out-band noise in dB;
        n_components = 2; %Number of components to estimate;
        yl=[-1 1]*10;
        xl = [30 50];
   case 'nsr_ecg_chi2_noise'  % Denoising of Normal ECG
        
        db = 'nsrdb'; %% Database to query
        recordID =   '16483'; %% Record ID to get
        InbandNoise = 4; % In-band noise in dB
        OutbandNoise = 7; %Out-band noise in dB;
        n_components = 3; %Number of components to estimate;
        yl=[-1 1]*10;
        xl = [30 40];
    case 'arrhythmia'
        
            db = 'mitdb'; %% Database to query
        recordID =   '201'; %% Record ID to get
        InbandNoise = -20; % In-band noise in dB
        OutbandNoise = -20; %Out-band noise in dB;
        n_components = 4; %Number of components to estimate;
        clim = [0 1];
        yl=[-1 1]*6;
        xl = [0 30] + 5*60;

end

if ~exist('physionetdb.m','file')
    if exist('wfdb/mcode','dir')
      addpath('wfdb/mcode/')  
    else
        fprintf(['This example attempts to retrieve ECG data from the Physionet database\n',...
             'using the WFDB Toolbox for Matlab and Octave, available at\n\n\t<a href="https://www.physionet.org/physiotools/matlab/wfdb-app-matlab/">https://www.physionet.org/physiotools/matlab/wfdb-app-matlab/</a>\n\n',...
             'Please install it and add it to the Matlab path.\n'])
     return
    end
end

%%% To retrieve other available records from the current database use this:
% records = physionetdb(db);


switch db
    case 'nsrdb'
        %Fs = 128;
        Fs = 125;
    case 'mitdb'
        Fs = 360;
    otherwise
        Fs = [];
end

switch example
    case {'nsr_ecg_noise','nsr_ecg_chi2_noise'}  % Denoising of Normal ECG
        ekg.nsamp = Fs*120; % Get two minutes of data
        ekg.recstart =Fs*5*60; %Starting 5 minutes into the recording
    
    case 'arrhythmia'
        ekg.nsamp = Fs*600; % Get ten minutes of data
        ekg.recstart =Fs*15*60; %Starting 15 minutes into the recording

end


[ecg,Fs,tm]=rdsamp(sprintf('%s/%s',db,recordID),1,ekg.nsamp+ekg.recstart,ekg.recstart,1);
  
%%% Standardize to unit variance
ecgz = zscore(ecg(:));

%%% Get signal power spectral density
sigSpectrum = abs(fft(ecgz).^2)./sqrt(ekg.nsamp);
whiteNoise = zscore(randn(size(ecgz)));

switch example
    case 'nsr_ecg_chi2_noise'
        whiteNoise = zscore(imag(hilbert(zscore(whiteNoise.^2))));
end

whiteNoiseFT = fft(whiteNoise);
%%% Create noise
NoiseInband  = zscore(real(ifft(whiteNoiseFT.*sigSpectrum./(1+sigSpectrum))));
NoiseOutband = zscore(real(ifft(whiteNoiseFT./(1+sigSpectrum))));
NoiseCombined = 10.^(InbandNoise/20)*NoiseInband + 10.^(OutbandNoise/20)*NoiseOutband;

ecgz_noise = ecgz+NoiseCombined;

%%% Get the window duration in samples
L = round(window_duration*Fs);

%%% Initialize the HOS object
clear hos;
hos(n_components) = hosobject(3);
hos.initialize(L,Fs,lowpass);

%%% Train on the input data through a maximum of 25 iterations
hos.get_block(ecgz_noise,25);
recovered_ecg = hos.xrec(ecgz_noise);
filtered_ecg = hos.xthresh(ecgz_noise);

%%% SNR improvement: Here measured as the difference between fisher Z transformed correlations before
%%% and after reconstruction, converted to dB (10*log10) units.
fisherzdB = @(x)5*log10(((1+x)./(1-x)));
cr = corr(ecgz,[ecgz_noise,recovered_ecg]);
snr_improvement = fisherzdB(cr)-fisherzdB(cr(1));
[mxcr,mxi] = max(cr(2:end));
[srt,srti] = sort(cr(2:end),'descend');

cols = 'rcmgy';
cols = cols(mod(srti-1,length(cols))+1);
%% Make a plot
t = (0:length(ecgz)-1)/Fs;

figure('units','normalized','position',[ 0.1333    0.1708    0.6547    0.6433])
subplot(4,1,1)
ylim(yl)
%xlim([30 50])
plh=[];
plh2=[];
switch example
    case {'nsr_ecg_noise','nsr_ecg_chi2_noise'}
        plh(end+1)= plot(t,ecgz_noise,'color',[1 1 1]*.5);
        title(sprintf('ECG + %idB in-band noise + %idB out-band noise',InbandNoise,OutbandNoise))
        legend({'Signal + Noise'})
         ylim(yl)
        pause(2)
end
hold on, 
plh(end+1)=plot(t,ecgz,'k','linewidth',2);
title('Original ECG')
legend(plh(end:-1:1),{'Original Signal','Signal + Noise'})
ylim(yl)


nplot = length(hos);
wb = fftshift(hos(1).freqindx.Bfreqs{1});

subplot(2*floor(nplot/2),ceil(nplot/(floor(nplot/2))),floor(nplot/2)*ceil(nplot/(floor(nplot/2)))+1)
    bc = hos(1).B./hos(1).D; BC = bc(hos(1).fullmap); %Partial bicoherence using the normalization of the full signal.
%BC = hos(1).bicoh; %Partial bicoherence on the residual signal with bias correction.
%Remap to the full range of frequencies.
imagesc(wb,wb,fftshift(abs(BC)));
colorbar
caxis(clim)
axis image xy
title('Initial bichorence')
    
pause(2)

subplot(4,1,1)
switch example
    case {'nsr_ecg_noise','nsr_ecg_chi2_noise'}
        title(sprintf('Recovered ECG:\nCorrelation %0.2f -> %0.2f  (%0.2fdB improvement) in cmp. %i',cr(1),cr(mxi+1),max(snr_improvement(mxi+1)),mxi))
     case 'arrhythmia'
           title('Recovered ECG')
end


lg = {'Signal + Noise','Original Signal'};
lg2={};
for k = 1:length(hos)-1;
    subplot(4,1,1)
    lg{end+1} = sprintf('Recovered component %i',k);
    plh(end+1)=plot(t,recovered_ecg(:,k),cols(k));
    if k==mxi
        set(plh(end),'LineWidth',2)
    end
    legend(plh(end:-1:1),lg(end:-1:1))
    
    subplot(4,1,2)
 
    lg2{end+1} = sprintf('Filtered for  component %i',k);
    plh2(end+1)=plot(t,filtered_ecg(:,k),cols(k));
     title('Filtered and thresholded')
    hold on
    legend(plh2(end:-1:1),lg2(end:-1:1))
    
    subplot(2*floor(nplot/2),ceil(nplot/(floor(nplot/2))),floor(nplot/2)*ceil(nplot/(floor(nplot/2)))+k+1)
    bc = hos(k+1).B./hos(1).D; BC = bc(hos(1).fullmap); %Partial bicoherence using the normalization of the full signal.
%    BC = hos(k+1).bicoh; %Partial bicoherence on the residual signal with bias correction.
    %Remap to the full range of frequencies.
    imagesc(wb,wb,fftshift(abs(BC)));
    colorbar
    caxis(clim)
       axis image xy
     title(sprintf('Residual bichorence after component %i',k))
    pause(1)
end

 subplot(4,1,1)
xlcurr = xlim;
nzoom = 25;
for k = 1:nzoom
  
   subplot(4,1,1)
   xlim(xl.*k/nzoom + xlcurr*(1-k/nzoom))
   ylim(yl)
   subplot(4,1,2)
   xlim(xl.*k/nzoom + xlcurr*(1-k/nzoom))
  % ylim(yl)
   drawnow
end

%legend(plh(end:-1:1),lg(end:-1:1))
