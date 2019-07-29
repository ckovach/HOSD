
function out = constructModel(Y,fs,varargin)

if isscalar(Y)
    N = Y;
    Y=[];
else
    N = length(Y);
end

i = 1;

out.Trange = [];
out.eventTimes = [];
out.eventTypes = [];
out.Tstart = 0;
out.modelType = 'binomial';
out.fs = fs;

while i < varargin
    
    switch(varargin{1})
        
        case {'Trange','modelType','Tstart','eventTimes','eventTypes','history'}
            out.(varargin{i}) = varargin{i+1};
            i=i+1;
        otherwise
            if ~ischar(varargin{i})
                error('Expecting a keyword string, found a %s object instead.',class(varargin{i}));
            end
            error('%s is an unrecognized keyword.',varargin{i});
    end
    i=i+1;
end