function [dprime,precision,recall,accuracy] = snr_prediction(skw,rate,thresh)

%  [dprime,precision,recall,accuracy] = snr_prediction(skw,rate,thresh)
%
% Predict SNR from skewness and rate based on some simplified
% assumptions. Specifically, it is assumed that the filter output contains a
% signal of fixed amplitude impulses in additive iid Gaussian noise. 
% 
% Input:
%
%   skw - skewness
%   rate - rate as feature probability per observation interval
%   thresh - detection threshold after standardization
% 
% [dprime,...] = snr_prediction(hos,x)
%
%   If the input is a hosobject, hos, and data segment, x, skewness, rate and
%   threshold are computed for the hosobject with the supplied data.
%
% Output:
%
%   dprime - SNR as the difference between the noise and signal means over
%            noise RMS amplitude.
%   precision - positive predictive value,  TPR/(TPR + FPR)
%   recall  -   sensitivity, TPR/(TPR + FNR)
%   accuracy -  accuracy computed as TPR/(TPR + FNR + FPR)
%
%   

% C. Kovach 2023


%The following estimates the separation of the signal and noise distributions based on
%filter skewness. The signal and noise means are assumed to average to 0,
%so that if the signal mean is mu, the noise mean is -rate/(1-rate)*mu.

%mu = nthroot((1-rate).^2./(rate.*(1-2*rate)).*skw,3); %Estimated standardized mean of the signal distribution

if isa(skw,'hosobject')
   hos = skw;
   x = rate;
   [xrec,xf,xthr] =  hos(1).xrec(x);
   skw = cumulant(xf,3);
   rate = nnz(xthr)/length(xthr);
   thresh = nthroot(hos(1).current_threshold(zscore(xf)),3);
  
else
    hos = false;
end

mu = nthroot((1-rate).^2.*skw,3)./sqrt(nthroot((rate.*(1-2*rate)).^2,3) - rate.*nthroot((1-rate).*skw.^2,3)); %Estimated standardized mean of the signal distribution

% var = (1+rate./(1-rate).*mu.^2);
%thresh = thresh*sqrt(var);


dprime = 1./(1-rate).*mu; %Estimated standardized separation.

% precision = @(thresh)(1-normcdf(thresh-dprime))*rate./((1-normcdf(thresh-dprime))*rate + (1-normcdf(thresh))*(1-rate));
% recall = @(thresh)1-normcdf(thresh-dprime);


% precision = (1-normcdf(thresh-dprime)).*rate./((1-normcdf(thresh-dprime)).*rate + (1-normcdf(thresh)).*(1-rate));
% recall = 1-normcdf(thresh-dprime);


TPR = (1-normcdf(thresh-mu)).*rate; %Predicted true positive rate
TNR = normcdf(thresh+rate./(1-rate).*mu).*(1-rate);     %Predicted true negative rate.
FPR = (1-normcdf(thresh+rate./(1-rate).*mu)).*(1-rate); %Predicted false positive rate
FNR = normcdf(thresh-mu).*rate;  %Predicted false negative rate

precision = TPR./(TPR + FPR);
recall = TPR./(TPR+FNR);

 
 accuracy = TPR./(TPR + FNR + FPR); % For the purpose of characterizing detection 
                                    % accuracy, this ignores true-negative rate, as
                                    % it would otherwise dominate the estimate
                                    % estimate and make converge to ~(1-rate) ~ 1 
                                    % at high thresholds. 
                                    
if isa(hos,'hosobject') && size(hos,2)>1
    [dpr,prec,rec,acc] = snr_prediction(hos(2:end),x-squeeze(xrec));
   
    dprime = [dprime,dpr];
    precision = [precision,prec];
    recall = [recall,rec];
    accuracy = [accuracy,acc];
    
end
    