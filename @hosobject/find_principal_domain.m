function [PD,Ws,Is,keep] = find_principal_domain(me,freqs,order,lowpass,highpass,mask,xlowpass,xhighpass,slowpass,shighpass,diagonal_slice,include_sig)

% Find the principal domain in a higher-order spectrum
%
% The principal domain represents a set of non-redundant regions in the
% higher order spectrum. Redundancy has two sources: (1)
% permutation of axes and (2) conjugate symmetry. By (1) we can include
% only regions such that 0<s1*W1<s2*W2<...<sk*Wk, where sj is a sign and by (2) we 
% can exclude regions with signature given by a sign reversal of an included region. 
% We can therefore proceed by including regions on the basis of (1) for each 
% non-redundant signature under sign reversal. In the case of the bispectrum, 
% there is only one non-redundant signature (+,+,-), hence (1) applies
% directly. In the case of the trispectrum, there are two non-redundant
% signatures:  (+,+,+,-) and (+,+,-,-), and so on. 
%
%
% INPUT: 
%
%    freqs: array of frequencies or 1xorder cell array of frequencies
%    order: order of the higher-order spectrum.
%
% OUTPUT:
%
%    PD: Logical array indicating the principal domain
%    Ws: Cell array of frequencies
%    Is: Cell array of indices.
%    keep: Logical vector over the full (order-1)-dimensional grid of the
%          first order-1 frequencies, true at the retained points.
%
% Candidate points. Every criterion below only removes points, and three of
% them keep only a thin slab of the grid: the diagonal slice (some pair of
% frequencies equal), the "OR" lowpass xlowpass (some |w_k| below it) and,
% for order > 3, the pair-sum lowpass slowpass (every frequency has a
% partner with |w_i + w_j| <= slowpass). Each is a union of a few linear
% constraints on the first order-1 frequencies, so instead of laying out
% the whole grid (dims^(order-1) points) and testing every point, the
% candidates are enumerated from the thinnest such family directly (see
% SLABINDICES: the last constrained axis is solved for, the others are
% enumerated), or taken from an explicit mask, whenever that is smaller
% than the grid. The criteria are then applied to the candidates exactly
% as before, so the result does not depend on which route was taken; the
% full grid remains the fallback.

% Copyright Christopher Kovach, University of Iowa, 2018, 2026

if nargin < 12 
    include_sig = [];
end
if nargin < 11 || isempty(diagonal_slice)
    diagonal_slice = false;
end
if nargin < 10 || isempty(shighpass)
    shighpass = min(highpass);
end
if nargin < 9 || isempty(slowpass)
    shighpass = max(lowpass);
end
    
if nargin < 7 || isempty(xlowpass)
    xlowpass = xlowpass*ones(1,order);
end
if nargin < 6 || isempty(mask)
    mask = true;
end

if nargin < 3 || isempty(order)
    if iscell(freqs)
        order = length(freqs);
    else 
        order = 3;
    end
end
if nargin < 4 || isempty(lowpass)
    lowpass = Inf*ones(1,order);
end
if nargin < 5 || isempty(highpass)
    highpass = zeros(1,order);
end

if isscalar(lowpass)
    lowpass = ones(1,order)*lowpass;
end

if isscalar(highpass)
    highpass= ones(1,order)*highpass;
end
%%    
if isnumeric(freqs)    
    freqs = repmat({freqs},1,order);
end

%%% For cross-polyspectra involving fewer signals than the specified order,
%%% assume that the last signal is repeated.
if length(freqs) < order
    freqs(end+1:order) = freqs(end);
end

dims = cellfun(@length,freqs(1:order-1));
ngrid = prod(dims);
dgtol = 1.5*me.sampling_rate/me.buffersize;
passtol = .25*me.sampling_rate/me.buffersize;
gridtol = 1e-6*min(cellfun(@(f)min(diff(sort(f(:)))),freqs(1:order-1)));

%%% Families of slabs whose union contains every point the corresponding
%%% criterion keeps. A slab is lo <= sum_q c(q)*w_q <= hi over the first
%%% order-1 frequencies; w_order = -sum(w_1..w_{order-1}) is implied.
families = {};
if diagonal_slice
    sl = {};
    for k = 1:order-1                          % |w_k - w_order| < dgtol
        c = ones(1,order-1); c(k) = 2;
        sl{end+1} = struct('c',c,'lo',-dgtol,'hi',dgtol); %#ok<*AGROW>
    end
    for k = 1:order-2                          % w_k == w_kk
        for kk = k+1:order-1
            c = zeros(1,order-1); c(k) = 1; c(kk) = -1;
            sl{end+1} = struct('c',c,'lo',0,'hi',0);
        end
    end
    families{end+1} = struct('kind','slabs','slabs',{sl});
end
if any(~isinf(xlowpass)) && all(~isinf(xlowpass(1:order)))
    sl = {};
    for k = 1:order-1                          % |w_k| < xlowpass(k)
        c = zeros(1,order-1); c(k) = 1;
        sl{end+1} = struct('c',c,'lo',-xlowpass(k),'hi',xlowpass(k));
    end
    c = ones(1,order-1);                       % |w_order| < xlowpass(order)
    sl{end+1} = struct('c',c,'lo',-xlowpass(order),'hi',xlowpass(order));
    families{end+1} = struct('kind','slabs','slabs',{sl});
end
if order > 3 && (max(shighpass) >= 0 || any(slowpass<lowpass))
    %%% Survivors of the pair-sum lowpass have a partner for w_1: itself,
    %%% one of w_2..w_{order-1}, or w_order (|w_1 + w_order| = |w_2 + ... + w_{order-1}|).
    S = max(slowpass)-passtol;
    sl = {};
    c = zeros(1,order-1); c(1) = 2;
    sl{end+1} = struct('c',c,'lo',-S,'hi',S);
    for j = 2:order-1
        c = zeros(1,order-1); c(1) = 1; c(j) = 1;
        sl{end+1} = struct('c',c,'lo',-S,'hi',S);
    end
    c = ones(1,order-1); c(1) = 0;
    sl{end+1} = struct('c',c,'lo',-S,'hi',S);
    families{end+1} = struct('kind','slabs','slabs',{sl});
end
if ~isscalar(mask)
    families{end+1} = struct('kind','mask');
end

%%% Pick the thinnest family; fall back to the full grid when none is
%%% appreciably smaller than it.
best = 0; bestcount = Inf;
for f = 1:length(families)
    switch families{f}.kind
        case 'mask'
            count = nnz(mask);
        case 'slabs'
            count = 0;
            for s = 1:length(families{f}.slabs)
                sb = families{f}.slabs{s};
                [~,cnt] = slabIndices(freqs,dims,sb.c,sb.lo,sb.hi,gridtol,true);
                count = count+cnt;
                if count >= bestcount
                    break
                end
            end
    end
    if count < bestcount
        best = f; bestcount = count;
    end
end

if best > 0 && bestcount < 0.75*ngrid
    switch families{best}.kind
        case 'mask'
            lin = find(mask(:));
        case 'slabs'
            lin = cell(1,length(families{best}.slabs));
            for s = 1:length(lin)
                sb = families{best}.slabs{s};
                lin{s} = slabIndices(freqs,dims,sb.c,sb.lo,sb.hi,gridtol,false);
            end
            lin = unique(vertcat(lin{:}));         % sorted: the grid's linear order
    end
    kpindx = false(ngrid,1);
    kpindx(lin) = true;
    Is = repmat({[]},1,order-1);
    [Is{:}] = ind2sub(dims,lin);
else
    lin = [];                                     % every grid point is a candidate
    freqsi = cellfun(@(x)1:length(x),freqs(1:end-1),'uniformoutput',false);
    Is = repmat({[]},1,order-1);
    [Is{:}] = ndgrid(freqsi{:});
    Is = cellfun(@(x)x(:),Is,'uniformoutput',false);
    kpindx = true(ngrid,1);
end
%%% The criteria below act on the candidates through keepc (one entry per
%%% candidate); the full-size keep is assembled from it once, at the end.
keepc = true(numel(Is{1}),1);


%%% The number of regions in the principal domain depends on the possible 
%%% signatures.
nsig = floor(order/2);
signatures = (-1).^(repmat(1:order,nsig,1)- 1>= order-repmat((1:nsig)',1,order));


Ws = {};
Wsum=0;

for k = 1:length(Is)
    
    Ws{k}(:,1) = freqs{k}(Is{k});
    Wsum  = Wsum+Ws{k};
    
end
Ws{order} = -Wsum;
% 
if diagonal_slice
    %%% Count the equalities among the frequencies of each candidate: with
    %%% the last frequency within dgtol, among the first order-1 exactly.
    dgs = 0;
    for k = 1:order-1
        dgs = dgs + (abs(Ws{k} - Ws{order}) < dgtol);
    end
    for k = 1:order-2
        for kk = k+1:order-1
            dgs = dgs + (Ws{k} == Ws{kk});
        end
    end
    keepc = keepc & dgs>=diagonal_slice;

    Wseq = 0;
    for k = 1:length(Ws)-1
        for kk = k+1:length(Ws)
           Wseq = Wseq + (abs(Ws{k} - Ws{kk}) < dgtol);
        end
    end

    keepc = keepc & Wseq>=diagonal_slice;      
% else
%     keep  = true;
end

if any(~isinf(xlowpass))
    keepxlp=false;
    for k = 1:length(Ws)
        keepxlp = keepxlp | abs(Ws{k})<xlowpass(k);
    end
     keepc = keepxlp & keepc;
%     keep = keepxlp & keep;
end

if any(xhighpass>0)
    keepxhp=false;
    for k = 1:length(Ws)
        keepxhp = keepxhp | abs(Ws{k})>=xhighpass(k);
    end
     keepc  = keepc & keepxhp;
%     keep  = keep & keepxhp;
end

for k = 1:length(Ws)
   
    if ~isinf(lowpass(k))
         keepc = keepc & abs(Ws{k})<lowpass(k)-passtol;
%         keep = keep & abs(Ws{k})<lowpass(k);
    end
    if highpass(k)>0
%         keep = keep & abs(Ws{k})>highpass(k);
         keepc = keepc & abs(Ws{k})>highpass(k)+passtol;
    end
end
if ~isscalar(mask)
    if isempty(lin)
        keepc = keepc & mask(:);
    else
        keepc = keepc & mask(lin);
    end
end

if ~isempty(include_sig)
    Wsig=abs(sum(sign([Ws{:}]),2));
    
    keepc = keepc & ismember(Wsig,include_sig);    
end

for k = 1:length(Ws)
     Ws{k} = Ws{k}(keepc);
%     Ws{k} = Ws{k}(keep);
    if k<=length(Is)
         Is{k} = Is{k}(keepc);
%         Is{k} = Is{k}(keep(kpindx));
    end
end




[~,srti] = sort(lowpass(1:end-1)); %This ensures the principal domain includes the complete range of frequencies
srti(end+1) = order;
PD = false;
for k = 1:nsig
    
    if ~isempty(include_sig) && ~ismember(abs(sum(signatures(k,:))),include_sig)
        continue
    end
    
    PD0 = Ws{srti(1)}>=0 & Ws{order}<=0; %First signature is always + and last always -.  ;
    for kk = 2:order
%           PD0 = PD0 & signatures(k,kk)*Ws{kk}>=signatures(k,kk-1)*Ws{kk-1};
        if signatures(k,srti(kk))==signatures(k,srti(kk-1)) 
           PD0 = PD0 & signatures(k,srti(kk))*Ws{srti(kk)}>=signatures(k,srti(kk-1))*Ws{srti(kk-1)};
        else            
           PD0 = PD0 & signatures(k,srti(kk))*Ws{srti(kk)}>=0;
        end
    end

    PD = PD | PD0;
  
end

if  order >3
    
    %%% Identify remaining regions that are not unique under conjugation.
    WW = [Ws{:}];
    pdi = find(PD);
%     [~,ismi] = ismember(sort(-WW(pdi,:),2),sort(WW(pdi,:),2),'rows');
    [~,ismi] = ismembertol(sort(-WW(pdi,:),2),sort(WW(pdi,:),2),'ByRows',true);
    [~,srti] = sortrows(WW(pdi,:));
    rnk(srti) = 1:length(srti);
    PD(pdi(ismi>0))= rnk(ismi(ismi>0)) <= rnk(ismi(ismi(ismi>0)));

    if  (max(shighpass) >= 0 || any(slowpass<lowpass))
        %%% Now we need to account for regions within the principal domain for which some
        %%% subset of the frequencies falls within the highcut range. This is
        %%% essentially the subset sum problem, which is NP complete. Here we will limit the search to
        %%% frequency pairs, thereby excluding any component of HOS that
        %%% weights by the power spectrum. The resultant windowing is not strictly
        %%% quasicumulant for orders > 5, but it removes a dominant
        %%% contribution from the power spectrum.
        WKP = [Ws{:}];
        discard = false(size(WKP,1),1);
        for k = 1:order-1
            discard = discard | any(abs(repmat(WKP(:,k),1,order-k)+WKP(:,k+1:order)) < max(shighpass)+passtol,2)   | all(abs(repmat(WKP(:,k),1,order)+WKP) > max(slowpass)-passtol,2)  ;       
        end

        PD(discard) = [];
        for k = 1:length(Ws)
            Ws{k} = Ws{k}(~discard);
            if k<=length(Is)
                Is{k} = Is{k}(~discard);
            end
        end
        keepc(keepc) = ~discard;
    end
end

%%% The retained points on the full grid
keep = kpindx;
keep(kpindx) = keepc;

end

function [lin,count] = slabIndices(freqs,dims,c,lo,hi,tol,countOnly)

% [lin,count] = slabIndices(freqs,dims,c,lo,hi,tol,countOnly)
%
% Linear indices into the grid of the K = numel(dims) frequency axes
% freqs{1..K} of the points (i_1,...,i_K) with
%
%     lo - tol <= sum_q c(q)*freqs{q}(i_q) <= hi + tol
%
% (a superset of the exact slab, for the caller to filter). The last axis
% with a nonzero weight is solved for: the other axes are enumerated, and
% for each of their combinations the admissible values on that axis form
% an interval, located in the sorted frequencies with DISCRETIZE and
% expanded into indices. Cost: prod(dims(others)) combinations plus the
% number of points found, never the whole grid. With countOnly the second
% output is returned without building lin.

K = numel(dims);
Q = find(c ~= 0);
qL = Q(end);
others = setdiff(1:K,qL);

if isempty(others)
    G = {};
    ncombo = 1;
    P = 0;
else
    ivec = arrayfun(@(d)(1:d)',dims(others),'uniformoutput',false);
    G = cell(1,numel(others));
    if isscalar(ivec)
        G{1} = ivec{1};
    else
        [G{:}] = ndgrid(ivec{:});
        G = cellfun(@(g)g(:),G,'uniformoutput',false);
    end
    ncombo = numel(G{1});
    P = zeros(ncombo,1);
    for m = 1:numel(others)
        q = others(m);
        if c(q) ~= 0
            f = freqs{q}(:);
            P = P + c(q)*f(G{m});
        end
    end
end

fL = freqs{qL}(:);
[fs,srt] = sort(fL);
b1 = (lo - P)/c(qL);
b2 = (hi - P)/c(qL);
lower = min(b1,b2) - tol;
upper = max(b1,b2) + tol;
edges = [-Inf; fs; Inf];
nlow = discretize(lower,edges) - 1;      % frequencies strictly below the interval
nhigh = discretize(upper,edges) - 1;     % frequencies at or below its top
counts = max(nhigh - nlow,0);
count = sum(counts);
if countOnly
    lin = [];
    return
end

combo = repelem((1:ncombo)',counts);
first = cumsum([0; counts(1:end-1)]);
pos = nlow(combo) + ((1:count)' - repelem(first,counts));
iL = srt(pos);
subs = cell(1,K);
for m = 1:numel(others)
    subs{others(m)} = G{m}(combo);
end
subs{qL} = iL;
lin = sub2ind(dims,subs{:});

end
