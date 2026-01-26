function out = ximp(me,xin,return_sparse,apply_window)


% out = ximp(me,xin,return_sparse,apply_window)
%
% Apply filtering and thresholding to input data with deflation and
% return an indicator at samples that are above threshold.
% This filters data using the current feature detection filter and applies  
% the current threshold to the result. The residual signal after removing the
% reconstructed component signal (see HOSOBJECT/RECONSTRUCT) is then passed
% to the next element in the hosobject array.  
% 
% Inputs: 
%   xin - input data as a column vector or a me.buffersize x N matrix
%   return_sparse - output is sparse if true (default).
%   apply_window - if the input is a me.buffersize x N matrix, then apply
%            the window specified in me.window to columns of xin. Default
%            is false.
% Output:
%   out - Same as me(1,:).xthresh>0
%
% See also HOSOBJECT/XTHRESH
%
% Copyright Christopher K. Kovach, University of Iowa 2018-2021

   if nargin < 2
       xin = me.dat;
   end
   if nargin < 3  || isempty(return_sparse)
       %Cannot do sparse output if it requires more than 2
       %dimensions
       return_sparse = size(xin,2) < 2 || length(me) < 2;
   end
   if nargin < 4 || isempty(apply_window)
      apply_window = false; 
   end
   out = me.xthresh(xin,return_sparse,apply_window)~=0;  
%            if length(me)>1
%                out = [out,me(2:end).ximp(in-me(1).xrec(in))];
%            end         

end
