 function [Gout,sgn] = partial_delay_filt(me,Xs,returnfull,use_sample_bispectrum,normalization)
            
% [Gout,sgn] = partial_delay_filt(me,Xs,returnfull,use_sample_bispectrum)
%
% Calculate the partial delay filter for input sample, optionally
% using the sample bispectrum.
%
% Inputs:
%   Xs - segmented input as buffersize x N matrix or cell array of length
%        1 to order containing buffersize x N data segments. 
%   returnfull - If false, return only values of the filter FFT within the passband,
%                otherwise return the full filter FFT. Default is true.
%   use_sample_bispectrum - Compute bicoherence from the sample(s) in Xs,
%               otherwise use pre-computed statistics. Default is false.
% Outputs:
%   Gout - FFT of partial delay filter(s) for each column of the input. If
%          returnfull is false, only FFT coefficients within the passband 
%          specified in the hosobject [me.lowpass, me.highpass] are returned (for the sake of efficiency). 
%
% Copyright Christopher K. Kovach, University of Iowa 2018-2021
    if nargin < 3 || isempty(returnfull)
        returnfull = true;
    end
    if nargin < 4 || isempty(use_sample_bispectrum)
        use_sample_bispectrum = false; % Uses precomputed statistics if false
    end
    if nargin< 5 || isempty(normalization)
        normalization = true;
    end
    if ~iscell(Xs)
        Xs  = {Xs};
    end

    %m = size(Xs{1},2);


    if isempty(Xs{1})
        return
    end

       %%% Adjust for lag
%             dt = atan2(imag(me.lag),real(me.lag))/(2*pi)*me.fftN;
%             delt = me.radw*dt;
%             delt(isnan(delt))=0;
%             delt = 0;
    for k = 1:length(Xs)
%                 FXs{k} = repmat(exp(-1i*delt),1,size(Xs{k},2)).*fft(Xs{k});
        FXs{k} = fft(Xs{k});
    end

    if length(FXs) == me.order
        FX = FXs{me.order};                
    elseif length(FXs)>1 && length(FXs)<me.order
        FXs = [FXs(ones(1,me.order-length(FXs)-1)),FXs,FXs(1)]; 
        FX = FXs{1};
    else
        FX = FXs{1};
    end

%             Xwin = Xin.*me.win;
%             FX = fft(Xwin);
   % FFX = 1;
%            FFXpart = ones([size(me.freqindx.Is,1),size(FX,2),me.order]);

    if ~use_sample_bispectrum || ~me(1).do_bsp_update
        FFX = conj(FX(me.freqindx.Is(:,me.order),:));
        FFXpart = {};
        FFXpart(1:me.order-1) = {FFX};
        FFXpart{me.order} = ones(size(FFX));

        for k = me.order-1:-1:1
            if k>length(FXs)
                FX = FXs{1};
            else
                FX = FXs{k};
            end
%                 FX(isnan(FX)) = 0;
            FXk = FX(me.freqindx.Is(:,k),:);
            for kk = setdiff(1:me.order,k) %%% Need multiple me.order symmetry regions for avg. partial cross-polyspectra
                FFXpart{kk} = FFXpart{kk}.*FXk;
            end
            FFX = FFX.*FXk;
        end
    else
        FFXpart = me.update_bispectrum(FX,true);
    end

    if normalization
        BC = me.B./(me.D+eps);
        bias = sqrt(me.BIASnum./(me.D.^2+eps));
        bias(isnan(bias))=0;
        BC = (abs(BC)-bias).*BC./(abs(BC)+eps);
        H = conj(BC./(me.D+eps));
    else
        BC = me.B; 
        H = conj(BC);
    end

    Gpart = 0;
%             Bcheck = 0;
%              HFcheck =  H(1:end-1).*FFX; %Sanity check to make sure the integration is correct.
    len = length(me.B)-1;
    for k = 1:length(FFXpart)

        if k<=length(me.Imats)
            I = me.Imats{k};
            Iconj = me.Iconjmats{k};
        else
            Iconj = integrator(me.freqindx.remap.*cast(me.freqindx.PDconj & me.freqindx.partialSymmetryRegions==k,class(me.freqindx.remap)),1,len,[0 len+1]);
            I = integrator(me.freqindx.remap.*cast(~me.freqindx.PDconj & me.freqindx.partialSymmetryRegions==k,class(me.freqindx.remap)),1,len,[0 len+1]);
            me.Imats{k}= I;
            me.Iconjmats{k} = Iconj;
        end
        HF = repmat(H(1:end-1),1,size(FFXpart{k},2)).*FFXpart{k};

        Gpart = Gpart + I*HF + Iconj*conj(HF);
%                 Bcheck = Bcheck + I*HFcheck + Iconj*conj(HFcheck);
    end
    if   me.check_sign 
        GFX = zeros(size(FX));
        GFX(me.keepfreqs{1},:) = Gpart(me.keepfreqs{1}(abs(me.freqs{1})<=me.lowpass(1)),:).*FX(me.keepfreqs{1},:);
        GX = real(ifft(GFX));
        sgn = sum(sign(GX).*(abs(GX)==max(abs(GX)))) ==1;
        Gpart=Gpart.*sgn;
    else
         sgn=1;
    end

    if returnfull
        Gout = zeros(size(Xs{1}));
        Gout(me.keepfreqs{1},:) = Gpart(me.keepfreqs{1}(abs(me.freqs{1})<=me.lowpass(1)),:);
    else
        Gout = Gpart(me.keepfreqs{1}(abs(me.freqs{1})<=me.lowpass(1)),:);
    end
end
