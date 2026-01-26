  function out = xprob(me,xin,return_sparse,apply_window)
 
% out = probh(me,xin,[return_sparse],[apply_window])
%
% This function returns the output of xthresh transformed into a feature
% probability based on the standard deviation of the filtered residual.
% 
% Inputs: 
%   xin - input data as a column vector or a me.buffersize x N matrix
%   return_sparse - output is sparse if true (default) and input is a column vector.
%   apply_window - if the input is a me.buffersize x N matrix, then apply
%            the window specified in me.window to columns of xin. Default
%            is false.
% Output:
%   out - feature probability vector, which is the output of HOSOBJECT/XTHRESH
%   normalized by residual std dev and passed through a standard gaussian CDF with v 
%
% See also HOSOBJECT/XTHRESH HOSOBJECT/IMP HOSOBJECT/APPLY_FILTER HOSOBJECT/RECONSTRUCT
% HOSOBJECT/CURRENT_THRESHOLD HOSOBJECT/FILTER_THRESHOLD 
%
% Copyright Christopher K. Kovach, University of Iowa 2018-2021

   if nargin < 2 || isempty(xin)
       xin = me.dat;
   end
   if nargin < 4 || isempty(apply_window)
      apply_window = false; 
   end
   if nargin < 3 || isempty(return_sparse)
       return_sparse = true;
   end
   %Cannot do sparse output if it requires more than 2
   %dimensions
   return_sparse = return_sparse & (size(xin,2) < 2 || length(me) < 2);

   xf = me(1).apply_filter(xin,apply_window,false);
%            if size(in,1)==me(1).buffersize %% Make sure the output is consistent if the input happens to be of buffersize length
%                xf = ifftshift(xf,1);
%            end


   if return_sparse
       xthr = sparse(double(me(1).filter_threshold(xf)));  
   else         
       xthr = double(me(1).filter_threshold(xf));       
   end
    mxthr = nanmean(xthr(xthr~=0));
    sdxthr = nanstd(xthr(xthr~=0));
    
    xfresid = me(1).xfilt(xin-me(1).xrec(xin));
    xfnorm = xthr./nanstd(xfresid);
    noise_likelihood = 1-normcdf(xfnorm(xthr~=0));
    if mod(me(1).order,2)==2 %Two sided for even orders
        noise_likelihood = min(noise_likelihood, 1-noise_likelihood)/2;
    end
    pevents = mean(xthr~=0); %Estimate of baseline probability
    % 
    % This assumes that the expectation of the positive distribution is the
    % observed value with variance equal to noise. A more accurate approach
    % will be to estimate the distribution under the positive case, but
    % this shortcut is used for now.
    out = xthr;
    out(xthr~=0) = 1 - noise_likelihood*(1-pevents)./(1+pevents*noise_likelihood); 
   if length(me)>1
       out= cat(sum(size(xin)>1)+1,out,me(2:end).xthresh(xin-me(1).xrec(xin,[],apply_window),return_sparse,apply_window));
   end


end