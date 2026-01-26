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

   xf = me(1).apply_filter(xin,apply_window,false);
%            if size(in,1)==me(1).buffersize %% Make sure the output is consistent if the input happens to be of buffersize length
%                xf = ifftshift(xf,1);
%            end


   if return_sparse
       out = sparse(double(me(1).filter_threshold(xf)));  
   else         
       out = double(me(1).filter_threshold(xf));       
   end


   if length(me)>1
       out= cat(sum(size(xin)>1)+1,out,me(2:end).xthresh(xin-me(1).xrec(xin,[],apply_window),return_sparse,apply_window));
   end


end