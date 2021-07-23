 function out = xrec(me,xin,thresh,apply_window,varargin)
 
% out = xrec(me,xin,[thresh],[apply_window])
%
% Reconstruct component signals within the input data with deflation. 
% This reconstructs the data using the current feature detection filter, 
% threshold, and feature waveform, then passes the residual to the next
% element in the hosobject array.  
% 
% Inputs: 
%   xin - input data as a column vector or a me.buffersize x N matrix
%   thresh - static cumulant threshold (see HOSOBJECT/CURRENT_THRESHOLD)
%   apply_window - if the input is a me.buffersize x N matrix, then apply
%            the window specified in me.window to columns of xin. Default
%            is false.
% Output:
%   out - Same as cat(d, me(1,1).reconstruct(xin,...), me(1,2:end).xrec(xresid))
%         where d is 2 if xin is a column vector and 3 if xin is a matrix.
%         and xresid = xin -  me(1,1).reconstruct(xin,...).
%
% See also HOSOBJECT/APPLY_FILTER HOSOBJECT/RECONSTRUCT
% HOSOBJECT/CURRENT_THRESHOLD HOSOBJECT/FILTER_THRESHOLD
%
% Copyright Christopher K. Kovach, University of Iowa 2018-2021

   if nargin < 2
       xin = me.dat;
   end
   if nargin < 3
       thresh = [];
   end
   if nargin < 4
       apply_window = false;
   end
   out = me(1).reconstruct(xin,thresh,apply_window,varargin{:}); 
   out(isnan(out)) = 0;
   if length(me)>1
       out =  cat(sum(size(xin)>1)+1,out,me(2:end).xrec(xin-out,thresh,apply_window,varargin{:}));
   elseif all(out(:)==0)
       sz = num2cell(size(xin));
       sz{sum(size(xin)>1)+1}=length(me);
       out(sz{:})=0;
   end
end 