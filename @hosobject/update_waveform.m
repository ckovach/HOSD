
function update_waveform(me,FXsh,initialize)

% update_waveform(me,FXsh,initialize)
% 
% Update the value of the feature waveform according to the learning rate
% and the shifted input segment.
% Input:
%   FXsh - FFT of an input segment shifted according to the detection filter.
%   initialize - true if learning is reinitialized, so that prior values
%          are discarded.
%
% Copyright Christopher K. Kovach, University of Iowa 2018-2021

   m = size(FXsh,2);
   if nargin < 3 || isempty(initialize)
       initialize = false;
   end

   if initialize
       lradj=1;
   else
       lradj = me.learningfunction(me.filter_adaptation_rate,m,1./me.filter_adaptation_rate);
   end

%             lradj = (1-(1-me.filter_adaptation_rate)^m);

%            me.wavefft = me.wavefft*(1-lradj) + mean(FXsh,2)*lradj; 
   if isempty(me.sampweight)
       smpw = ones(size(FXsh,2),1)/size(FXsh,2);
   else
       smpw = me.sampweight;
   end
   me.wavefft = me.wavefft*(1-lradj) + (FXsh*smpw)*lradj; 

end