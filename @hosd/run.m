function [XF,xthr,betas] = run(me,x,start_at_component,varargin)

%Run HOSD fitting until stopping criteria are met.
%
% [XF,xthr,betas] = run(me,x) also returns each component's filter output
% (XF), thresholded filter output (xthr) and reconstruction scaling (betas),
% as HOSOBJECT/XREC lists them: one column per component for a record (a
% column vector x), one page per component, along the third dimension, when
% x holds windows as the columns of a matrix.

if me.blank 
    compi = 1;
elseif nargin < 3 || isempty(start_at_component)
    compi = length(me.hos)+1;
else
    compi = start_at_component;
end
if ~islogical(me.apply_ICA) & size(me.apply_ICA,1) == size(x,2)
    me.ICA_unmixing = me.apply_ICA;
    me.apply_ICA = true;
elseif ~islogical(me.apply_ICA)
    error('apply_ICA must be either a logical or an unmixing matrix matching the dimension of the input')
end
if me.apply_ICA && size(me.ICA_unmixing,1) ~= size(x,2)
    fi = which('fastica');
    if ~isempty(fi)
        if me.order == 3
            g = 'skew';
        else
            g = 'pow3';
        end
        [~,~,UM] = fastica(x(~any(isnan(x),2),:)','g',g);
        me.ICA_unmixing = UM';
    else
        UM = pica(x,[],me.order);
        me.ICA_unmixing = UM;
    end
elseif ~me.apply_ICA
    me.ICA_unmixing =1;
end

start_at_component=compi;
xresid = x*me.ICA_unmixing;



cml = Inf*ones(1,me.lookahead+1);

% Parameter template for the component(s) this call appends. Normally
% every component copies hosica(1), so the whole array shares one
% parameter set. A component_template (see hosd/addcomp) overrides that
% for this call, allowing a component with different parameters -- e.g.
% a different order -- to be appended; the residual chaining below (and
% in hosobject/xrec, /xdetect) is per element, so mixed chains are fine.
% Consumed here so it applies only to the components added by this run.
if ~isempty(me.component_template)
    hos0 = me.component_template;
    me.component_template = [];
else
    hos0 = me.hosica(1);
end

plh = me.liveplot;

while max(cml) > me.threshold && compi-start_at_component < me.maxcomps

    if isa(hos0,'mvhosd') || me.multivariate
        hos = mvhosd(hos0);

        if compi==1
            me.hosf = hosobject(hos0);
        end
    else
        hos = hosobject(hos0);
        if compi == 1
            me.hosf = hos;
        else
            me.hosf(compi) = hos;
        end
    end
    if compi > 1
        me.hosica(compi) = hos;
    else
        me.hosica = hos;
    end
    if me.multivariate
        [~,~,plh] = me.hosica(compi).get_block(xresid,me.niter,plh,[],compi);
    else
        [~,~,~,~,~,plh] = me.hosica(compi).get_block(xresid,me.niter,plh,[],compi);
    end
    me.liveplot = plh;
    [xr,xf,xt,beta]  = me.hosica(compi).xrec(xresid);
    %  xf = me.hos(compi).xfilt(xresid);
    if nargout > 0
        % The filter output has the input's dimensions, and the components
        % follow them: one column per component for a record, one page per
        % component for windows (the columns of a matrix).
        sub = repmat({':'},1,max(1,sum(size(xf)>1)));
        XF(sub{:},compi) = xf; %#ok<AGROW>
        xthr(sub{:},compi) = xt; %#ok<AGROW>
        betas(:,compi) = beta; %#ok<AGROW>
    end
    if me.multivariate
        me.hosf(compi) = hosobject(hos0);
        me.hosf(compi).get_input(xf);
    end

    % The standardized cumulant of the component's filter output over every
    % sample of it. For windows that is all the windows together: cumulant
    % works along the first dimension, so on the matrix it would return one
    % value per window.
    cml = [cumulant(xf(:),me.hosica(compi).order),cml];


    xresid = xresid-squeeze(xr);

    cml(end) = [];


    fprintf('\nComp. %i cumulant %0.2f',compi,cml(1));

    compi = compi+1;
end
me.blank = false;