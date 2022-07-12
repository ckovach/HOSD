function FFXpart = update_bispectrum(me,FXs,initialize)
            
    % Right now this updates HOS in chunks with temporal decay weighting
    % applied serially across chunks. That is, HOS is computed as a weighted
    % average of the simple within-chunk average and the prior value,
    % where weighting is determined by the integral of the serial weight
    % for the chunk.

    if ~iscell(FXs)
        FXs  = {FXs};
    end

    m = size(FXs{1},2);

    if nargin < 3 || isempty(initialize)
        initialize = false;
    end


    if isempty(FXs{1})
        return
    end

       %%% Adjust for lag
    dt = atan2(imag(me.lag),real(me.lag))/(2*pi)*me.fftN;
    delt = me.radw*dt;
    delt(isnan(delt))=0;
    for k = 1:length(FXs)
        FXs{k} = repmat(exp(-1i*delt),1,size(FXs{k},2)).*FXs{k};
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
    FFX = zeros(size(me.freqindx.Is,1),size(FX,2),size(FX,3));
    FXk = FFX;
    FFX(:) = conj(FX(me.freqindx.Is(:,me.order),:));

    FFXpart = {};
    FFXpart(1:me.order-1) = {FFX};
    FFXpart{me.order} = ones(size(FFX));

    for k = me.order-1:-1:1
        if k>length(FXs)
            FX = FXs{1};
        else
            FX = FXs{k};
        end
       
        FXk(:) = FX(me.freqindx.Is(:,k),:);
        for kk = setdiff(1:me.order,k) %%% Need multiple me.order symmetry regions for avg. partial cross-polyspectra
            FFXpart{kk}(:) = FFXpart{kk}.*FXk;
        end
        FFX = FFX.*FXk;
    end

    if isempty(me.sampweight)
%                 wgt = ones(size(FFX,2),1)/size(FFX,2);
         if me.integrated_magnitude_normalization
            ims = nansum(abs(FFX))';
            wgt = 1./(ims+nanmean(ims(:)).*1e-6);
            wgt(isnan(wgt))=0;
         else
             wgt = ~any(isnan(FFX))./sum(~any(isnan(FFX)),2);
         end
    end
    for kk =1:me.order
            FFXpart{kk}(isnan(FFXpart{kk})) = 0;
    end
    FFX(isnan(FFX)) = 0;

%             BX = mean(FFX,2);
    BX = sum(FFX.*wgt,2);
    XPSD = sum(wgt.*abs(FX).^2,2);
%             BXpart = mean(FFXpart,2);

    BX(end+1,:) = nan;
    XPSD(end+1,:) =nan;
%             BXpart(end+1,:) = 0;
    BXpart = {};
    for kk = 1:me.order
%                BXpart{kk} = mean(FFXpart{kk},2); 
       BXpart{kk} = sum(FFXpart{kk}.*wgt,2); 
       BXpart{kk}(end+1,:) = 0;
    end

    if initialize
        lradj = 1;
        fflr = 1;
        lrbias = 1;
        lr=1;
        me.sumlr = 1;
        me.sumlr2 = 1/size(FX,2);

    else
        [lradj,lr] = me.learningfunction(me.hos_learning_rate,m,me.hos_burnin);
         fflr = me.learningfunction(me.filter_adaptation_rate,m,me.filter_burnin);
%               fflr = (1-(1-me.filter_adaptation_rate)^m);
        %%% Adjust the learning rate according to the number of samples
        %%% in Xin, giving the sample weight the same total as if it
        %%% had been added serially.
         asympedf = 2./lr-1; %Asymptotic EDF

        lrbias = 1./asympedf*(1-(1-lr).^(2*m)); % "Learning rate" for the sum of squared weights in the bias term
        me.sumlr = me.sumlr*(1-lradj) + lradj;
        me.sumlr2 = me.sumlr2*(1-lr).^(2*m) + lrbias;
    end            
%     BX(isnan(BX))=0;
    me.B = (1-lradj)*me.B + lradj*BX;
%             me.Bpart = (1-fflr)*me.Bpart + fflr*BXpart;
     for kk = 1:me.order
        me.Bpart{kk} = (1-fflr)*me.Bpart{kk} + fflr*BXpart{kk};
     end
%     XPSD(isnan(XPSD))=0;
    me.PSD = (1-lradj)*me.PSD + lradj*XPSD;
%            me.sumlr2 = (me.sumlr2-1./asympedf)*(1-me.current_learning_rate).^(2*m) + lrbias;

    switch me.normalization
        case 'awplv'
%                     NX = mean(abs(FFX),2);
            NX = sum(abs(FFX).*abs(wgt),2);
            NX(end+1,1) = 0;
%                     XbiasNum = sum(abs(FFX).^2,2)./m^2;
            XbiasNum = sum(abs(FFX).^2.*abs(wgt).^2,2);
            XbiasNum(end+1,1) = 0;
            XbiasNum(isnan(XbiasNum))=0;
            NX(isnan(NX))=0;
            me.BIASnum = me.BIASnum.*(1-lr).^(2*m) + lrbias*XbiasNum;
            me.D = (1-lradj)*me.D + lradj*NX+eps;

        case {'bicoh','bicoherence'}
            %This doesn't seem to have strictly correct symmetry
            XBCpart = sum(wgt.*abs(FFXpart{1}).^2,2);
            XBCpart(end+1,1) = 0;
            XBCpart(isnan(XBCpart))=0;
            me.BCpart = (1-lradj)*me.BCpart + lradj*XBCpart;                    
            me.D = sqrt(me.BCpart.*me.PSD(me.freqindx.Is(:,1)))+eps;
    end

end