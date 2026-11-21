function [XF,xthr,betas] = addcomp(me,x,varargin)

%%% Estimate the next component
%
% addcomp(me,x)              Append a component using the current
%                            settings (a copy of hosica(1)'s parameters).
% addcomp(me,x,template)     Append a component whose parameters come
%                            from the hosobject (or mvhosd) TEMPLATE
%                            instead of hosica(1). The appended component
%                            may differ in order or band settings from
%                            the preceding ones -- e.g. a 4th-order
%                            component fit to the residual of a 3rd-order
%                            one -- because the hosobject array chains
%                            residually per element.
%
% In every form the new component is fit to the residual left by the
% components already in the array (see the xrec subtraction below).

tmpl = [];
if ~isempty(varargin) && ...
        (isa(varargin{1},'hosobject') || isa(varargin{1},'mvhosd'))
    tmpl = varargin{1};
    varargin(1) = [];
end

if me.blank
    me.component_template = tmpl;
    if nargout > 0
        [XF,xthr,betas] = me.run(x,varargin{:});
    else
        me.run(x,varargin{:});
    end
else
    % The sum of the components' reconstructions. xrec lists the components
    % along the first dimension after the input's own (see HOSOBJECT/XREC):
    % along 2 for a record (a column vector), along 3 for windows (the
    % columns of a matrix), whose own second dimension is the windows. A
    % multivariate fit lists them along 2, with the channels along 3.
    xr = me.hos.xrec(x);
    if isa(me.hosica,'mvhosd')
        xr = squeeze(sum(xr,2));
    else
        xr = sum(xr,sum(size(x)>1)+1);
    end
    % Set the template only once the residual is in hand: if xrec throws,
    % a template parked on the handle would otherwise survive and be
    % silently consumed by the next, unrelated run().
    me.component_template = tmpl;
    if nargout > 0
        [XF,xthr,betas] = me.run(x-xr,varargin{:});
    else
        me.run(x-xr,varargin{:});

    end
end