 function [lradj,lr] = learningfunction(me,learningrate,m,burnin)

% [lradj,lr] = learningfunction(me,learningrate,m,burnin)
%  
% Calculate the sample weighting based on the learning rate, current window 
% number and burn-in parameter and current sample size. When m > 1,
% the updated estimate is a weighted average of the prior estimate and the
% sample estimate, where the latter is weighted according to the integral
% of the learning function if each of the m samples had been added serially. 
%
% Inputs:
%   learningrate - Baseline learning rate between 0 and 1, where 0 
%   m   -   Current input sample size. Learning rate is adjusted to account
%           for the sample size.
%   burnin - Burnin period. Learning rate is initially 1./(me.window_number+1)
%            and relaxes to the baseline learning rate with time constant 
%            given by burnin. If burnin is infinite  then sample weighting is 
%            equivalent to a simple uniformly weighted average over samples.
% Outputs:
%   lradj - learning rate adjusted for both window number (according to burnin) 
%           and sample size.
%   lr - learning rate adjusted for window number but not sample size.
%
% See also HOSOBJECT/DO_UPDATES
%
% Copyright Christopher K. Kovach, University of Iowa 2018-2021
   

    if nargin < 3 || isempty(m)
        m = 1;
    end
%     if nargin < 4 || isempty(burnin)
%         burnin = me.burnin;
%     end
    lr = exp(-me.window_number./burnin)./(me.window_number+1) + (1-exp(-me.window_number./burnin))*learningrate;
    lradj = (1-(1-lr)^m);

end