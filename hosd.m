
classdef hosd < handle

% This class implements HOSD with a cumulant threshold for choosing the number of components. 

properties
    hos
    order
    maxcomps = 10;
    threshold = .05;
    multivariate = false;
    niter = 25;
    liveplot = true;
    lookahead = 0;
    hosf = [];

end

methods

    function me = hosd(varargin)
        
        if nargin > 0
        if isa(varargin{1},'mvhosd')
           
            me.hos = mvhosd(varargin{1});
            me.multivariate = true;

        elseif isa(varargin{1},'hosobject')
            
            me.hos = hosobject(varargin{1});
            me.multivariate = false;
        end
        end
        
        props = fieldnames(me);
        rm = [];
        for k = 1:length(props)
            idx = find(strcmpi(varargin,props{k}));
            if ~isempty(idx)
                rm = [rm,idx,idx+1];    
                me.(props{k}) = varargin{idx+1};
            end
        end
        varargin(rm) = [];

        if isempty(me.hos)
            if me.multivariate
                me.hos = mvhosd(varargin{:});
            else
                me.hos = hosobject(varargin{:});
            end
        end
        me.order = unique([me.hos.order]);

    end


    function run(me,x,varargin)
        
        compi = 1;
        xresid = x;

        cml = Inf*ones(1,me.lookahead+1);
    
        hos0 = me.hos(1);
  
        plh = me.liveplot;

        while max(cml) > me.threshold && compi < me.maxcomps
            
            if isa(hos0,'mvhosd') || me.multivariate
               hos = mvhosd(hos0);
               
               if compi==1
                   me.hosf = hosobject(hos0);
               end
            else
                hos = hosobject(hos0);

            end
            if compi > 1
                me.hos(compi) = hos;
            else
                me.hos = hos;
            end
            [~,~,plh] = me.hos(compi).get_block(xresid,me.niter,plh,[],compi);
            xf = me.hos(compi).xfilt(xresid);
            
            if me.multivariate
                me.hosf(compi) = hosobject(hos0);
                me.hosf(compi).get_input(xf);
            end

            cml = [cumulant(xf,me.hos(compi).order),cml];
            
            xr  = me.hos(compi).xrec(xresid);

            xresid = xresid-squeeze(xr);
            
            cml(end) = [];
                

            fprintf('\nComp. %i cumulant %0.2f',compi,cml(1));

            compi = compi+1;
        end
      
    end
end

end