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
   if length(me)>1 && ismatrix(xin) && size(xin,1) > me(1).fftN
       %%% A record deflated for the components after this one: the
       %%% reconstruction filters it on the way, and that filter output is
       %%% this component's. RECONSTRUCT returns it with NaN read as 0, so
       %%% NaN goes back wherever the FIR reaches a non-finite input sample
       %%% (the output is kept ceil(fftN/2) samples after its input; see
       %%% APPLY_FILTER); there the output is NaN or +/-Inf and never finite,
       %%% and RECONSTRUCT keeps the Inf. Windows (fftN or fewer rows) take
       %%% the plain path: RECONSTRUCT fftshifts the filter output of fftN-row
       %%% windows.
       [xr,out] = me(1).xrec(xin,[],apply_window);
       xnf = any(~isfinite(xin),2);
       if any(xnf)
           out(movmax(xnf,[floor(me(1).fftN/2)-1, ceil(me(1).fftN/2)]) & ~isinf(out)) = NaN;
       end
       xr(isnan(xr)) = 0; % as XREC with one output (all NaN for a component with no supra-threshold sample)
   else
       [out,~] = me(1).apply_filter(xin,apply_window,false);  
       if length(me)>1
           xr = me(1).xrec(xin,[],apply_window);
       end
   end

   if length(me)>1
       out = cat(sum(size(xin)>1)+1,out,me(2:end).xfilt(xin-xr,apply_window));
   end
end