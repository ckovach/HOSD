
classdef ppmodel 
%
%  Regressor object for point-process models. 
% 
%  out=ppmodel(Response,RegIn,varargin)
%
%  Response is a struct with two fields:
%       .dat : Response data
%       .type: data type : 'binary' for binary data,'count' for count
%                          data,'times' for event times.
%       .fs : sampling rate.
% 
%  RegIn is a struct describing the independent variables.
%       .autodep: How to model the conditional dependence on the process history
%                 using a laguerre polynomial expansion. It should be a struct
%                 array with the following fields:
%                    .order: order of the expansion.
%                    .tau:   time constant.
%                 An empty array indicates no modeling of conditional dependence;
%                 with no dependence, statistical results are valid only if the
%                 process is Poisson. Default is
%                 .autodep = struct('order',{0 1 2},'tau',{0 .1 .5});

% C. Kovach 2019

    properties

       response;
       labels;
       modeltype= 'binomial';
    end

    properties( Dependent = true)
       designMtx 

    end

    methods

        function me = ppmodel(me,varargin)
            if nargin > 1
              
                me.response = varargin{1};
            end
                
        end
        
        function out=get.designMtx(me)
            % Get the design matrix
        end
        
        function addregressor(me,regIn)
           
            
            
        end
        
    end
end

