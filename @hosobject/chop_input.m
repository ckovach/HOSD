 function [Xchop,T] = chop_input(me,xin,apply_window,delay)
 
% [Xchop,T] = chop_input(me,xin,[apply_window],[delay])
%
% Chop the input into windows of me.buffersize duration with overlap according to
% me.poverlap.
%
% Inputs:
%   xin - Input data as a column vector.
%   apply_window - If true, apply window specified in me.window after
%          segmentation.
%   delay - Adjust the timing of each window according to the value(s) in delay. 
%
% Outputs:
%   Xchop - segmented data.
%   T - segmentation matrix for the input.
%
% Copyright Christopher K. Kovach, University of Iowa 2018-2021

    if nargin < 3 || isempty(apply_window)
        apply_window = true;
    end
    if nargin < 4 || isempty(delay)
        delay=0;
    end
    nxin = length(xin);
    stepn = round(me(1).poverlap*me(1).bufferN);
    nget = nxin - me(1).bufferN+1;
    tindx = (0:me(1).bufferN-1)';
    wint = (0:stepn:nget-1)+delay;

    T = repmat(tindx,1,length(wint))+repmat(wint,length(tindx),1)+1;
    T(T>length(xin))=length(xin);
    T(T<1)=length(xin);
    Xchop = xin(T);
    if apply_window
       Xchop = fftshift(repmat(me(1).win,1,size(T,2)).*Xchop,1);
    end

end