
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

    %%%%% functions defined in separate files
    [XF,xthr,betas] = addcomp(me,x,varargin)        
    [XF,xthr,betas] = run(me,x,start_at_component,varargin)

end


end