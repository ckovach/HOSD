function out = xfilt(me,xin,apply_window)

% out = xfilt(me,in,apply_window)
%
% Apply filtering to input data with deflation. This filters the data using
% the current feature detection filter and passes the residual to the next element 
% in the hosobject array.  
% 
% Inputs: 
%   xin - input data as a column vector or a me.buffersize x N matrix
%   apply_window - if the input is a me.buffersize x N matrix, then apply
%            the window specified in me.window to columns of xin. Default
%            is false.
% Output:
%   out - Same as [ me(1,1).apply_filter(xin,...), me(1,2:end).xfilt(xresid)]
%         where xresid = xin -  me(1,1).reconstruct(xin,...).
%
% See also HOSOBJECT/APPLY_FILTER HOSOBJECT/RECONSTRUCT
%
% Copyright Christopher K. Kovach, University of Iowa 2018-2021


   if nargin < 2
       xin = me.dat;
   end
   if nargin < 3 || isempty(apply_window)
      apply_window = false; 
   end
    [out,~] = me(1).apply_filter(xin,apply_window,false);  

   if length(me)>1
       out = cat(sum(size(xin)>1)+1,out,me(2:end).xfilt(xin-me(1).xrec(xin,[],apply_window),apply_window));
   end
end