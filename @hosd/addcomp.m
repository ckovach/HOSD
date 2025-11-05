function [XF,xthr,betas] = addcomp(me,x,varargin)

%%% Estimate the next component

if me.blank
    if nargout > 1
        [XF,xthr,betas] = me.run(x,varargin{:});
    else
        me.run(x,varargin{:});
    end
else
    xr = squeeze(sum(me.hos.xrec(x),2));
    if nargout > 1
        [XF,xthr,betas] = me.run(x-xr,varargin{:});
    else
        me.run(x-xr,varargin{:});

    end
end