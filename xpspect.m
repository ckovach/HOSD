function [out,pspindices] = xpspect(Xs,fs,bws,order,varargin)

%[psp,pspindices] = xpspect({X1,...,Xk},fs,bws,order,[options])
%
% General function for computing polyspectra and cross polyspectra of the 
% given  order using the window-overlap method.
%
% Input arguments:
%       Xs:  1 x order cell array with either column vectors or matrices
%          containing time x interval samples
%       fs:  sampling rate.
%       bws: vector of bandwidths to use for each element of Xs
%   order:  Order of the polyspectrum.
%
%
%Options: Options may either be specified as pairs of keywords and values, 
%         i.e. xpspect(...,'keyword',value,...) or as one or more structs
%         with  struct.keyword = value.
%   Avalailable options are:      
%
%   lowpass: limit the range of frequencies for each axis to values less-
%           than or equal to this. This may be specified as a scalar or as
%           a vector of order-1 length, which applies a separate limit for
%           each dimension.
%   highpass: limit to frequencies above this value, etc. 
%   maxfreq: Sum of frequencies across all dimensions must be less-than-or
%            -equal to this value.
% Output arguments:
%       psp: Struct with the following fields:
%          .pspect: Normalized (default) or unnormalized polyspectrum as an (order-1)-dimensional
%                   matrix, depending on the "return_normalized" option.
%          .fs:     Frequuency labels for the dimensions of pspect
%          .options Options struct.
%          .segment Struct used in segmenting the data. 
%       pspindices: Struct containing indices into the original data with the
%            following fields:
%           .findex: index of the frequency for each unqiue term in the
%                    estimate.
%           .conjugate: Terms for which the complex conjugate is taken.
%           .reconmat: reconstruct into the same shape and size as psp.pspect. 
%           .lin     Vector of indices into pspect to create a vector of unique coefficients
%                    (inverse of reconmat).
%
%
%
% See also PSPECT2

% C. Kovach 2022

if nargin < 4
    order = 3;
elseif ischar(order)
    varargin = [{order},varargin];
    order = 3;
end

%%% OPTIONS & DEFAULT VALUES%%%
options.lowpass= Inf;
options.maxfreq= Inf;
options.highpass= 0;
options.normalization = 'awplv';
options.full_range = false; % Add negative frequencies if they are not already included
options.min_range = false; % Return only one symmetry region
options.symmetrize = false;
options.round_freq = true; % Round to the nearest frequency band if necessary.
options.tolerance = []; % Rounding tolerance (defaults to min(diff(f))).
options.stats = false; %Compute parametric statistics on the mean.
%options.real_signal=true;
options.window = 'sasaki';
options.upsampleFx = ones(1,order);
options.segment = [];
options.povlp = .5;
options.return_normalized = true;% The pspect field contains normalized and bias corrected HOS if true.
                                 % Unnormalized HOS is then (psp.pspect+psp.bias)*psp.normalization
% options.window_spacing = 'max'; %Space windows according to the longest segment.
% options.window_spacing = 'min'; %Space windows according to the shortest segment
options.window_spacing = 'mean'; %Space according to mean segment duration.

if isstruct(bws)
    options.segment = bws;
else
    options.segment = struct('fs',fs);
end

optfld = fieldnames(options);
i = 1;
while i <length(varargin)
    if isstruct(varargin{i})
        fldn = fieldnames(varargin{i});
        for k = 1:length(fldn)
            if ~ismember(fldn{k},optfld)
                error('Unrecognized option, %s',fldn{k})
            end
            options.(fldn{k}) = varargin{i}.(fldn{k});
        end
        i = i-1;
    elseif ismember(varargin{i},optfld)
        options.(varargin{i})=varargin{i+1};
    else
        error('Unrecognized option, %s',varargin{i});
    end
    i = i+2;
end
segin = options.segment;

if length(bws) < order
    bws(end+1:order) = bws(end);
end

if length(options.upsampleFx) < order
    options.upsampleFx(end+1:order) = options.upsampleFx(end);
end

if isnumeric(Xs)
    X = Xs;
    Xs = {};
    if size(X,3)>1
        X = permute(X,[1 3 2]);
    end
    for k = 1:size(X,2)
        Xs{k} = squeeze(X(:,k,:)); %#ok<*AGROW>
    end
    clear X
end

if length(Xs)<order
    Xs(end+1:order) =Xs(end);
end


%%% If input is not already chopped, then chop
if ~isfield(options.segment,'windur')
    options.segment.windur = ceil(fs./bws)/fs;
elseif length(options.segment.windur) < order
    options.segment.windur(end+1:order) = options.segment.windur(end);
end
lens = cellfun(@(x)size(x,1),Xs);

if ~isfield(options.segment,'wint')
    switch options.window_spacing
        case 'max'
            dwin = max(options.segment.windur); %Space according to the longest window. 
        case 'min'
            dwin = min(options.segment.windur); %Space according to the shortest window. 
        case 'mean'
            dwin = mean(options.segment.windur);%Space according to average window duration.
    end
    wint = dwin/2:dwin*options.povlp:min(lens)/fs-dwin/2;
    options.segment.wint = wint;
end
if ~isfield(options.segment,'Tranges')
        options.segment.Tranges = [-.5 .5]'*options.segment.windur;
end


discard = false;
for k = 1:length(Xs)
    if size(Xs{k},2)>1
        winN = ceil(size(Xs{k},1)*(options.upsampleFx(k)+1));
        X = Xs{k};
    else
        winN = ceil(options.segment.windur(k)*(options.upsampleFx(k)+1)*fs);
    

        segment = options.segment;
        segment.Trange = segment.Tranges(:,k)';
        segment.maxN = length(Xs{k});
        T = chopper(segment);
        
        X = Xs{k}(T).*window(options.window,size(T,1));
    end
    X(end+1:winN,:) = 0;
    FXs{k} = fft(X).';
    discard = discard | any(isnan(X));
    ws{k} = ifftshift((0:winN-1)-floor(winN/2))/winN*fs;
    FXs{k} = FXs{k}(:,ws{k}>=0);
    ws{k}(ws{k}<0) = [];
end

if any(discard)
    fprintf('\n%i segments with nan''s discarded (%0.2f%%)',sum(discard),mean(discard)*100)
    for k = 1:length(FXs)
        FXs{k}(discard,:) = [];
    end
end


[out,pspindices] = pspect2(FXs,ws,order,options);
out.segment = options.segment;
out.options.segment = segin;

if options.return_normalized
    out.pspect = (abs(out.pspect)./out.normalization - out.bias).*out.pspect./(abs(out.pspect)+eps);
end

