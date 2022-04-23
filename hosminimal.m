
classdef hosminimal < hosobject
   
    % A version of hosobject that stores minimal data needed to do signal reconstruction
    % and to recreate a full hos object minus the actual HOS. 
    % No HOS data are stored, only filters and waveforms.
    % Usage: 
    %   To create a hos object
    %        obj = hosminmal(hos)
    %   where hos is a hosobject.
    % 
    % To recover a hos object with the same parameters (but HOS statistics unset) 
    % use hosobject: hos = hosobject(obj);
    %
    % SEE HOSOBJECT
    
   
  
    
    
  
    properties (Dependent = true)
%         B
%         D
%         Imat
%         Iconjmat
    end
    properties
        subspace_dim=0;
        projection=1;
    end
    methods
       
        function me = hosminimal(order,varargin)
            
            if nargin ==0
                return
            end
            if isa(order,mfilename) || isa(order,'hosminimal') || isa(order,'struct')|| isa(order,'hosobject')
               obj = order;
%                fns = [{'BIASnum'};setdiff(properties(obj),{'BIAS','Bfull','H','bicoh','current_threshold','sampling_rate','freqindx','buffersize','filterftlag','fullmap','partialbicoh','filterfft','filterfun','bicohreduced'})];
              fns = setdiff(fieldnames(obj),{'freqs','BIAS','Bfull','H','bicoh','current_threshold','freqindx','filterftlag','fullmap','partialbicoh','bicohreduced','B','D','Imats','Iconjmats','BIASnum','freqindex','Bpart'});
             
               metac = metaclass(me);
               props = metac.PropertyList;

               getprops = strcmp({props.SetAccess},'public');
               fns = intersect(fns,{props(getprops).Name}); %This ensures that only fields with public set access are set to avoid unexpected behavior.
               
               if isa(order,'struct')
                   fnsunset = setdiff({'fftN','bufferN','sampling_rate','lowpass','freqs','freqindx'},fns);
                   for k = 1:length(fnsunset)
                       for kk = 1:length(obj)
                           obj(kk).(fnsunset{k})=me(1).(fnsunset{k});
                       end
                   end
               end               
               if length(obj)==1
                   obj(2:length(me)) = obj;
               elseif length(obj)>length(me)
                   me(length(obj)) = hosobject;
               end
               
               
               me(1).order = obj(1).order; 
               if obj(1).check_sign
                    me(1).check_sign = mod(obj(1).order,2)~=0;
               end
%                me.initialize(obj(1).bufferN,obj(1).sampling_rate,obj(1).lowpass,obj(1).freqs,obj(1).freqindx,varargin{:})
               
               me(1).do_indexing_update = false;
               priority = intersect({'order','buffersize','pad','lag','sampling_rate','keepfreqs'},fns);% These fields should be set first
               for k = 1:length(priority)
                   me(1).(priority{k}) = obj(1).(priority{k});
               end
               fns = setdiff(fns,priority);
               for k = 1:length(fns)  
                   if isprop(me(1),fns{k})
                        try
                       me(1).(fns{k}) = obj(1).(fns{k});
%                        if me(1).lag~=obj(1).lag
%                            keyboard
%                        end
                        catch
                        end
                   end
               end
               me(1).do_indexing_update = false;
               if length(me)>1
                   me(2:end) = hosminimal(obj(2:end));
               end
               return
            elseif  me(1).check_sign
               me(1).check_sign = mod(order,2)~=0;
            end
            
            if nargin < 1 
                return
            elseif nargin == 1
                me.order = order;
                return
            else
                me.order = order;
            end
            
%             me.initialize(varargin{:});
        end
        
        function reset(me)
          return
        end
        function initialize(me,varargin)
            return
        end
%         function out = get.B(me)     
%             out = [];
%         end
%         function out = get.D(me)     
%             out = [];
%         end
%         function set.B(me,~)     
%             warning('This object does not store HOS data')
%         end
%         function set.D(me,~)     
%             warning('This object does not store HOS data')
%         end
        
        function get_input(me)
            warning('This object does not implement HOSD')
        end
        function get_block(me)
            warning('This object does not implement HOSD')
        end
    end
    
end
