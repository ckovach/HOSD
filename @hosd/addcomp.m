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
    if nargout > 1
        [XF,xthr,betas] = me.run(x,varargin{:});
    else
        me.run(x,varargin{:});
    end
else
    xr = squeeze(sum(me.hos.xrec(x),2));
    % Set the template only once the residual is in hand: if xrec throws,
    % a template parked on the handle would otherwise survive and be
    % silently consumed by the next, unrelated run().
    me.component_template = tmpl;
    if nargout > 1
        [XF,xthr,betas] = me.run(x-xr,varargin{:});
    else
        me.run(x-xr,varargin{:});

    end
end