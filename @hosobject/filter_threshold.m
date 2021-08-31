function [Xthresh,trialthresh] = filter_threshold(me,Xfilt,thresh,use_adaptive_threshold)
            
% [Xthresh,trialthresh] = filter_threshold(me,Xfilt,[thresh],[use_adaptive_threshold])
% 
% Apply a threshold to the output of the feature detection filter. Thresholding is
% set such that the residual value of the static cumulant of order me.order
% is not greater than a specified value (default is 0).
%
% Inputs:
%     Xfilt - Data segment(s) to which the detection filter has been
%            applied
%     thresh - Value of the static cumulant according to which the
%              threshold will be set (default is 0).
%     use_adaptive_threshold - If false the cumulant is calculated
%            directly from sample, otherwise it is based on a reference 
%            distribution (see HOSOBJECT/CURRENT_THRESHOLD).
%
%Outputs:
%     Xthresh - Thresholded data.
%     thrialthresh - threshold value for the current input. Thresholding is
%               given by Xfilt^me.order >= trialthresh.
%
%
% See also HOSOBJECT/APPLY_FILTER HOSOBJECT/CURRENT_THRESHOLD
%
% Copyright Christopher K. Kovach, University of Iowa 2018-2021


    % Apply a moment-based threshold
    if nargin < 3 || isempty(thresh)
        thresh= me.thresh;
    end
    if nargin < 4 || isempty(use_adaptive_threshold)
        use_adaptive_threshold = me.use_adaptive_threshold;
    end

%                   zsc = @(x)(x-nanmean(x))./nanstd(x);
%                   Xcent = zsc(Xfilt);

    if me.use_adaptive_threshold
         Xcent = (Xfilt-me.running_mean)./sqrt(me.running_var);
    else
        Xcent = (Xfilt - nanmean(Xfilt(:)))./nanstd(Xfilt(:));
    end
    %Xcent = Xfilt;
    if isempty(me.threshold_order)
        me.threshold_order = me.order;
    end
    Xmom = Xcent.^me.threshold_order;

%             if size(Xfilt,1) == me.bufferN && size(Xfilt,2)==1 && use_adaptive_threshold
    if use_adaptive_threshold && ~isempty(me.CDFbuffer)%|| size(Xfilt,1) == me.bufferN && size(Xfilt,2)==1
         trialthresh = me.current_threshold([],thresh);
%                  Xcs = [];
    else
        trialthresh = me.current_threshold(Xcent,thresh);

    end

    if mod(me.threshold_order,2) < mod(me.order,2) %If using kurtosis for threshold setting with odd order HOSD, account for sign.
        Xmom(Xcent<0)=0;
    end

    switch me.threshold_type
        case 'hard'
            THR = Xmom>=trialthresh;%repmat(trialthresh,size(Xmom,1),1);
        case 'soft'
            sigm = @(x)1./(1+exp(-x));
            THR = sigm(me.threshtemp*(Xmom-repmat(trialthresh,size(Xmom,1),1)));
        otherwise
            error('Unrecognized threshold type')
    end
    Xthresh = Xfilt.*THR;
    Xthresh(isnan(Xthresh))=0;
end