
function out=createRegressors(Response,RegIn,varargin)

%
%  out=createRegressors(Response,RegIn,varargin)
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


