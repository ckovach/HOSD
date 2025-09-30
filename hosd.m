
classdef hosd < handle

% This class implements HOSD with a cumulant threshold for choosing the number of components. 

properties
    order
    maxcomps = 1;
    threshold = .05;
    multivariate = false;
    niter = 25;
    liveplot = true;
    lookahead = 0;
    hosf = [];
    blank = true;
    apply_ICA = false; %Use ICA at outset
    ICA_unmixing = [];
    hosica
end

properties (Dependent = true)
    hos
end

methods

    function me = hosd(varargin)
        
        if nargin > 0
            if isa(varargin{1},'hosd')
                obj = varargin{1};
                me.hosica = hosobject(obj.hos);
                props = setdiff(fieldnames(obj),'hos');
                for k = 1:length(props)
                    me.(props{k}) = obj.(props{k});
                end

            else
                
                if isa(varargin{1},'mvhosd')
                   
                    me.hosica = mvhosd(varargin{1});
                    me.multivariate = true;
        
                elseif isa(varargin{1},'hosobject')
                    
                    me.hosica = hosobject(varargin{1});
                    me.multivariate = false;
                  elseif isa(varargin{1},'hosd')
                    me.hosica = hosobject(varargin{1});
        
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
        
                if isempty(me.hosica)
                    if me.multivariate
                        me.hosica = mvhosd(varargin{:});
                    else
                        me.hosica = hosobject(varargin{:});
                    end
                end
                me.order = unique([me.hosica.order]);
            end
            if ~islogical(me.apply_ICA) && max(size(me.apply_ICA))>1
                me.ICA_unmixing = me.apply_ICA;
                me.apply_ICA = true;
            end
        end

    end

    function [XF,xthr,betas] = addcomp(me,x,varargin)
       
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
            
    end

    function [XF,xthr,betas] = run(me,x,start_at_component,varargin)
        
        
        if me.blank 
            compi = 1;
        elseif nargin < 3 || isempty(compi)
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
    
        hos0 = me.hosica(1);
  
        plh = me.liveplot;

        while max(cml) > me.threshold && compi-start_at_component < me.maxcomps
            
            if isa(hos0,'mvhosd') || me.multivariate
               hos = mvhosd(hos0);
               
               if compi==1
                   me.hosf = hosobject(hos0);
               end
            else
                hos = hosobject(hos0);

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
            [xr,xf,xthr,beta]  = me.hosica(compi).xrec(xresid);
          %  xf = me.hos(compi).xfilt(xresid);
            if nargout > 1
                XF(:,compi) = xf;
                xthr(:,compi) = xthr;
                betas(:,compi) = beta;
            end
            if me.multivariate
                me.hosf(compi) = hosobject(hos0);
                me.hosf(compi).get_input(xf);
            end

            cml = [cumulant(xf,me.hosica(compi).order),cml];
            
          
            xresid = xresid-squeeze(xr);
            
            cml(end) = [];
                

            fprintf('\nComp. %i cumulant %0.2f',compi,cml(1));

            compi = compi+1;
        end
        me.blank = false;
    end
%%%%%%
    function out = get.hos(me)
        if ~isempty(me.apply_ICA) && me.apply_ICA && ~me.blank && isa(me.hosica,'mvhosd') 
            for k = 1:length(me.hosica)
                out(k) = mvhosd(me.hosica(k));
                out(k).filterfun = permute(squeeze(out(k).filterfun)*me.ICA_unmixing',[1 3 2]);
                out(k).feature = permute(squeeze(out(k).feature)*pinv(me.ICA_unmixing),[1 3 2]);
            end
        elseif isa(me.hosica,'mvhosd'   )         
            out = mvhosd(me.hosica);
        else
            out = hosobject(me.hosica);
        end
    
    end

end

end