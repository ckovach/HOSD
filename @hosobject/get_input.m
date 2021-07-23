function out = get_input(me,xin,apply_window,use_shifted,initialize)

% X = get_input(me,xin,[apply_window],[use_shifted],[initialize])
%
% Implements online HOSD by updating the current values of the HOS, filter
% and feature estimates given an input data segment.
%   1. Updates the current estimates for the bispectrum, bicoherence normalization
%      provided me.update_bispectrum is true.
%   2. Aligns the input segment(s) according to extrema in with the current detection
%   filter.
%   3. Computes the partial delay filter(s) for the aligned segment. 
%   4. Updates the feature detection filter, and feature waveform, provided
%      me.update_filter and me.update_waveform are true, respectively.
%   5. Implements realtime deflationary HOSD by calling
%            me(1,2:end).get_input(xresid,...) 
%      where 
%            xresid = xin - me(1,1).reconstruct(xin)
% Updated values are computed as a weighted average of prior values and
% sample estimates according to the learning_rate parameters, me.hos_learning_rate
% and me.filter_adaptation_rate, me.burnin and me.window_number. 
%
% Note that realtime HOSD is non-iterative (that is, there is only 1 iteration
% in the realignment of each segment, while convergence relies on the
% serial alignemnt of input segments). For iterative offline HOSD, use
% me.get_block
%
% Inputs:
%   xin - input data as a column vector or me.buffersize x N matrix.
%   apply_window - Apply the window specified in me.window to input data
%          segments. Default is true.
%   use_shifted - Shift the input according to extrema in the detection
%          filter. Default is true.
%   initialize - Re-initialize estimates, discard all prior value.
%
% Output:
%   X - Segmented input data. 
%
%
% See also HOSOBJECT/DO_UPDATES HOSOBJECT/GET_BLOCK HOSOBJECT/RECONSTRUCT
% HOSOBJECT/CHOP_INPUT
%
% Copyright Christopher K. Kovach, University of Iowa 2018-2021
    nxin = numel(xin);
      if nargin < 5 || isempty(initialize)
         initialize=false;
      end
      if nargin < 4 || isempty(use_shifted)
         use_shifted=true;
     end
     if nargin < 3 || isempty(apply_window)
         apply_window=true;
     end
    if nxin >= me(1).bufferN
        me(1).bufferPos = 0; % Discard the buffer
        if  size(xin,1) ~= me(1).bufferN 
%                     stepn = round(me(1).poverlap*me(1).bufferN);
%                     nget = nxin - me(1).bufferN+1;
%                     tindx = (0:me(1).bufferN-1)';
%                     wint = (1:stepn:nget);
% 
%                     T = repmat(tindx,1,length(wint))+repmat(wint,length(tindx),1);
%                     Xchop = xin(T);
            [Xchop,T] = me(1).chop_input(xin,false);

            snip = xin(T(end)+1:numel(xin));
            if ~isempty(snip)
                me(1).write_buffer(snip);
            end
        else
            Xchop = xin;
        end
        me(1).do_updates(Xchop,apply_window,use_shifted,initialize)

    else
        me(1).write_buffer(xin);
        return
    end    
    out = Xchop;
    if length(me)>1
       xrec = me(1).reconstruct(xin,[],[],true);
       xrec(isnan(xrec)&~isnan(xin))=0;

       out = [out,me(2:end).get_input(xin-xrec,apply_window,use_shifted,initialize)];
    end



end