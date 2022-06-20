 function out =current_threshold(me,Xcent,thresh)
       
% out =current_threshold(me,[Xcent],[thresh])
%        
% Compute the current thresholding for the output of the feature detection filter. 
% The output is thresholded such that the value of the static cumulant of
% order me.order, computed on subthreshold samples, is not greater than the specified value
% (default 0). The threshold can be computed for a given sample or using
% the reference distribution contained in me.CDFbuffer. 
%
% Inputs:
%   Xcent - input data. If not given or empty, the threshold is determined
%      from me.CDFbuffer.
%   thresh - static cumulant thresholod. By default this is 0.
%
% Outputs:
%   out - Threshold value for Xcent.^me.order. To obtain a threshold value
%         for Xcent (without the exponent), use thresh = nthroot(out,me.order).
%         
%   
% Copyright Christopher K. Kovach, University of Iowa 2018-2021

if nargin < 2 || isempty(Xcent) %||true % Find threshold based on empirical CDF
    Xcent =  me.CDFbuffer(:,1);
    Xcent2 = me.CDFbuffer(:,2);
    XcentK = me.CDFbuffer(:,me(1).threshold_order);
else %Find threshold based on input moment
    Xcent = sort(Xcent(:));


    Xcent2 = Xcent.^2;
    XcentK = Xcent.^me.threshold_order;
end
if nargin<3 || isempty(thresh)
    thresh = me.thresh;
end
if isempty(me.threshold_order)
    me.threshold_order=me.order;
end
 if me.threshold_order == 3
    % For the bispectrum compute normalized skewness
%                 keepsamples = ones(size(Xcent));

%                   srt = sort(Xcent(:));
%                  outlier_threshold = 5;
     if me.outlier_threshold~=0
        keepsamples = ~isnan(iterz(Xcent,me.outlier_threshold,-1)); % Suppress extreme negative outliers           
     else
         keepsamples = ~isnan(Xcent);
     end
    Xcent(isnan(Xcent))=0;
    Xcent2(isnan(Xcent2))=0;
    XcentK(isnan(XcentK))=0;

    m1 = cumsum(Xcent.*keepsamples)./cumsum(keepsamples); % cumulative mean on sorted peaks
    m2 = cumsum(Xcent2.*keepsamples)./cumsum(keepsamples); % cumulative 2nd moment
    m3 = cumsum(XcentK.*keepsamples)./cumsum(keepsamples); % cumulative 3rd moment
    %  Third cumulant
    c3 = m3 - 3*m2.*m1 + 2*m1.^3; % Third cumulant on sorted peaks

    keepsrt = Xcent>0 & c3>  thresh;
    out = sum ((diff(keepsrt)>0).*XcentK(2:end,:));
    if ~any(out)
        out=Inf;
    end
else
    % For now, apply a simple threshold on the standardized moment for
    % orders > 3. This should be improved to use the proper
    % cumulant.
%                 Xmom = Xcent.^me.order;
    Xmom = XcentK;
    [Xsrt,srti] = sort(Xmom);
%                 Xpow = Xcent(srti).^2;
    Xpow = Xcent2(srti,:);
    Xmean = Xcent(srti,:);
%                 Xcent(isnan(Xcent))=0;
%                 Xcent2(isnan(Xcent2))=0;
%                 XcentK(isnan(XcentK))=0;
    Xpow(isnan(Xpow)) = 0;
    Xmean(isnan(Xmean))=0;

    keepsamples = ones(size(Xsrt));
    Mcs = cumsum(Xsrt.*keepsamples)./cumsum(keepsamples);
    Powcs = cumsum(Xpow.*keepsamples)./cumsum(keepsamples);
    Mncs = cumsum(Xmean.*keepsamples)./cumsum(keepsamples);
    if mod(me.threshold_order,2)==0
        %%% Correction for power spectral component with even
        %%% orders
        Xbaseline = (me.threshold_order-1)*(Powcs -   Mncs).^(me.threshold_order./2) +  Mncs.^me.threshold_order + thresh;
    else
        Xbaseline =thresh;
    end
    Mstd = (Mcs- Xbaseline)./Powcs.^(me.threshold_order/2) ;
    Xthr = Mstd>thresh;
    Xthr = cumsum(diff([zeros(2,size(Xthr,2));Xthr])>0,'reverse')==0; % Use the last threshold crossing if there are multiple
    if all(Xthr)
        out = Inf;
    else
        out = nansum(Xsrt.*(diff(Xthr)>0),1);
    end
 end
end