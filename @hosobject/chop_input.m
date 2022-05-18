 function [Xchop,T,segment] = chop_input(me,xin,apply_window,delay,segment)
 
% [Xchop,T,segment] = chop_input(me,xin,[apply_window],[delay],[segment])
%
% Chop the input into windows of me.buffersize duration with overlap according to
% me.poverlap.
%
% Inputs:
%   xin - Input data as a column vector.
%   apply_window - If true, apply window specified in me.window after
%          segmentation. Note that output is fftshifted along the 1st
%          dimension if apply_window is true.
%   delay - Adjust the timing of each window according to the value(s) in delay. 
%   segment - A struct with fields wint 
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
    if nargin < 5 || isempty(segment) || isstruct(segment) && isempty(segment.wint)
        nxin = length(xin);
        stepn = round(me(1).poverlap*me(1).bufferN);
        nget = nxin - me(1).bufferN+1;
        wint = (0:stepn:nget-1)+delay;
        segment.wint = wint;
    end
    if islogical(segment)
        if segment
            segment = me(1).segment;
        else
            error('Segment must be a structure or a logical value and true')
        end
    end
        
    if ~isfield(segment,'fs')
        segment.fs = 1;
    end
    
    if ~isfield(segment,'Trange')
        segment.Trange = [0 (me(1).bufferN-1)/segment.fs];
    end
    
%     tindx = (0:me(1).bufferN-1)';
    tindx =(segment.Trange(1):1/segment.fs:segment.Trange(2))';
    
    if apply_window && (~isfield(segment,'window') || isempty(segment.window))
        
        if isa(me(1).window,'function_handle')
            segment.window = me(1).window(length(tindx));
        else
            segment.window = window(me(1).window,length(tindx));
        end
        
    else
        segment.window = rectwin(length(tindx));
    end
    
    wint = segment.wint;
  
    T = round((repmat(tindx,1,length(wint))+repmat(wint,length(tindx),1))*segment.fs)+1;
    T(T>length(xin))=length(xin);
    T(T<1)=length(xin);
    Xchop = xin(T);
    if apply_window
       Xchop = fftshift(repmat(me(1).win,1,size(T,2)).*Xchop,1);
    end

end