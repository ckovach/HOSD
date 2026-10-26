  function out = xthresh(me,xin,return_sparse,apply_window)
 
% out = xthresh(me,xin,[return_sparse],[apply_window])
%
% Apply filtering and thresholding to input data with deflation. 
% This filters data using the current feature detection filter and applies  
% the current threshold to the result. The residual signal after removing the
% reconstructed component signal (see HOSOBJECT/RECONSTRUCT) is then passed
% to the next element in the hosobject array.  
% 
% Inputs: 
%   xin - input data as a column vector or a me.buffersize x N matrix
%   return_sparse - output is sparse if true (default) and input is a column vector.
%   apply_window - if the input is a me.buffersize x N matrix, then apply
%            the window specified in me.window to columns of xin. Default
%            is false.
% Output:
%   out - Same as cat(d, me(1,1).filter_threshold(xin,...), me(1,2:end).xthresh(xresid))
%         where d is 2 if xin is a column vector and 3 if xin is a matrix.
%         and xresid = xin -  me(1,1).reconstruct(xin,...).
%
% See also HOSOBJECT/APPLY_FILTER HOSOBJECT/RECONSTRUCT
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

   if length(me)>1 && ismatrix(xin) && size(xin,1) > me(1).fftN
       %%% A record deflated for the components after this one: the
       %%% reconstruction filters and thresholds it on the way (see XFILT).
       [xr,~,xthr] = me(1).xrec(xin,[],apply_window);
       xr(isnan(xr)) = 0; % as XREC with one output (all NaN for a component with no supra-threshold sample)
   else
       xthr = me(1).filter_threshold(me(1).apply_filter(xin,apply_window,false));
       if length(me)>1
           xr = me(1).xrec(xin,[],apply_window);
       end
   end

   if return_sparse
       out = sparse(double(xthr));  
   else         
       out = double(xthr);       
   end


   if length(me)>1
       out= cat(sum(size(xin)>1)+1,out,me(2:end).xthresh(xin-xr,return_sparse,apply_window));
   end


end