function out = freq2index(me,freqsin,order,lowpass,highpass,keepfreqs,condense,frequency_spacing,mask,xlowpass,xhighpass,slowpass,shighpass,diagonal_slice,include_sig)

% out = freq2index(freqs,order,lowpass,highpass)
%
% Function to generate indexing for polyspectra of a given order.
%
% INPUT:
%
%   freqs - vector or 1 x order cell array of vectors containing frequencies.
%   lowpass - lowpass frequency
%   highpass - highpass frequency
%
% OUTPUT:
%
%   out - struct with fields:
%       .Is : indices into the signal spectrum to produce a non-redundant
%             set of coefficients in the HOS of given order.
%       .remap : mapping from the non-redundant vector to the full HOS
%                array.
%       .principal_domain: Indicates the principal domain. 
%       .PDconj: Indicates regions related to the principal domain by
%                complex conjugate.

% Copyright Christopher Kovach, University of Iowa 2018.

if nargin < 15 || isempty(include_sig)
    include_sig = [];
end
if nargin < 14 || isempty(diagonal_slice)
    diagonal_slice = false;
end
if nargin < 13 || isempty(shighpass)
    shighpass = [];
    
end
if nargin < 12 || isempty(slowpass)
    slowpass = [];
end
if nargin < 11 || isempty(xhighpass)
    xhighpass = 0;
end
if nargin < 10 || isempty(xlowpass)
    xlowpass = Inf;
end

if nargin < 9 || isempty(mask)
    mask = true;
end
if nargin < 8 || isempty(frequency_spacing)
    frequency_spacing = 'linear';
end
if nargin < 3 || isempty(order)
    if iscell(freqsin)
        order = length(freqsin);
    else 
        order = 3;
    end
end
if nargin < 4 || isempty(lowpass)
    lowpass = .5;
end
if nargin < 5 || isempty(highpass)
    highpass = 0;
end

if isscalar(lowpass)
    lowpass = ones(1,order)*lowpass;
end
if isscalar(xlowpass)
    xlowpass = ones(1,order)*xlowpass;
end

if isscalar(highpass)
    highpass= ones(1,order)*highpass;
end
if isscalar(xhighpass)
    xhighpass = ones(1,order)*xhighpass;
end

%%
if isnumeric(freqsin)
    freqsin = repmat({freqsin},1,order);
    if nargin < 7 || isempty(condense)
        condense = true;
    end
elseif nargin < 7 || isempty(condense)
    condense = false; % If the indexing is the same for the different frequencies, then we will give the mapping into the principal domain
end
if nargin < 6 || isempty(keepfreqs)
   keepfreqs = cellfun(@(x)1:length(x),freqsin,'uniformoutput',false);
elseif isnumeric(keepfreqs)
    keepfreqs = repmat({keepfreqs},1,order);
elseif length(keepfreqs)<order
    keepfreqs(end+1:order) = keepfreqs(end);
end



%%% For cross-polyspectra involving fewer signals than the specified order,
%%% assume that the last signal is repeated.
if length(freqsin) < order
    freqsin(end+1:order) = freqsin(end);
end


keeplp = arrayfun(@(fr,lp)abs(fr{1})<=lp,freqsin(1:end-1),lowpass(1:end-1),'uniformoutput',false);

keeplp2 = cellfun(@(kpfr,kplp)kpfr(kplp),keepfreqs(1:end-1),keeplp,'uniformoutput',false);
dims = cellfun(@(x)sum(x),keeplp);
keepregion = false(dims);
keepregion(keeplp2{:})=true;

%%% Check if negative frequencies are explicitly represented
%%% Will index into the negative frequencies if so.
freqs = cellfun(@(fr,kp)fr(kp),freqsin,keepfreqs,'uniformoutput',false);

nneg = sum(freqs{end}<0);
npos = sum(freqs{end}>0);
two_sided = nneg>npos/2; %%% Assume that negative frequencies are just padding if there aren't at least as many negative frequencies as half the number of positive frequencies.

if ~isscalar(mask)
    mask = mask(keepregion);
end

[PD,Ws,Is,keep] = find_principal_domain(me,freqs,order,lowpass,highpass,mask,xlowpass,xhighpass,slowpass,shighpass,diagonal_slice,include_sig);
Fsum = Ws{end};

%%% Efficiently map the nearest elements of Fsum to elements of frinds{end} with
%%% the sort function. Frequencies in the range freq(1)-df(1)/2 will be
%%% mapped to freq(1) and freq(end)+df(end)/2 will be mapped to freq(end).
n = numel(Fsum);
[frsrt,frsrti] = sort(freqs{order});
switch frequency_spacing
    case 'linear'
        dfr = diff(frsrt)/2;
        frcent = dfr([1:end,end]) + frsrt;
        frcent = [-dfr(1) + frsrt(1),frcent];
    case {'log','logarithmic'}
        dlfr = diff(log(frsrt))/2;
        lfrcent = dlfr([1:end,end]) + log(frsrt);
        lfrcent = [-dlfr(1) + log(frsrt(1)),lfrcent];
        frcent = real(exp(lfrcent));
        frcent(isnan(frcent))=0;
end

%   frsrti(end+1) = length(freqs{order})+1;
[srt,srti] = sort([-Fsum(:);frcent(:)]);
E = [zeros(n,1);ones(length(frcent),1)];
E = E(srti); 
IND = zeros(size(E));
IND(srti) = cumsum(E);
 IND(n+1:end) = [];
%  IND(IND==0) = length(freqs{end});
IND(IND==0) = 1;
IND(IND>length(frsrti)) = length(frsrti);
IND =  frsrti(reshape(IND,size(Fsum)));

Is{order}=IND(:) ;

% keep = IND>0 & IND<=length(freqs{end});

freqindex = cellfun(@(x)[find(x),0],keepfreqs,'uniformoutput',false);
Is = cellfun(@(x,fri)fri(x(:))',Is,freqindex,'uniformoutput',false);
Is = [Is{:}];
%Isreduced = Is(keep,:);
% W = [];
% for k = 1:order
%    W(k,:) = Ws{k};%(keep);
% end
W = [Ws{:}]';

% Use the most economical integer type for indexing
intTypes = {'uint8','uint16','uint32','uint64'};
maxInts = cellfun(@(x)double(intmax(x)),intTypes);
useInt = intTypes{logical(diff([false,size(Is,1)<maxInts]))};

subremap = zeros(size(keep),useInt);
subremap(keep) = find(keep);

tol = min(abs(diff([freqs{:}])));
tol = 2.^floor(log2(tol)-5); %Use power of 2 to avoid precision error
if condense %for auto-spectra we only need the principal domain. This is not so for cross-spectra
    %%% Canonical form of every kept tuple: magnitudes sorted ascending (a
    %%% negative value sorts after a positive one of equal magnitude),
    %%% rounded to tol, with the original signs.
    [Wsrt,wsrti] = sort(abs(W)+tol/32*(W<0)); %For frequencies of equal magnitude and different sign, force consistent sorting of negative and positive values
    Wsrt = round(abs(Wsrt)./tol)*tol;
    %[Wsrt,wsrti] = sort(round(abs(W)./tol)*tol);
    wsrti = wsrti + order*repmat(0:size(W,2)-1,order,1);
    Wsrt = Wsrt.*sign(W(wsrti));

    %%% Each kept point maps to the principal-domain tuple with the same
    %%% canonical form (ism) or the negated one (ismconj). The row-wise
    %%% lookups are done on one integer key per tuple when the keys fit in
    %%% a uint64, which is the same comparison as ismember(...,'rows') on
    %%% the rounded values (ismember returns the lowest matching index in
    %%% both forms), without sorting multi-column matrices.
    [ism,ismi] = rowsMember(Wsrt,Wsrt(:,PD),tol);
    [ismconj,ismiconj] = rowsMember(-Wsrt,Wsrt(:,PD),tol);
   ismiconj(1) = 0;
   ismconj(1)=false;
   ismiconj(ism & ismconj)=0; %Ignore if both are in the principal domain.
   ismconj(ism & ismconj) = 0;
%   IsreducedPD = Isreduced(PD(keep),:);
   IsPD = Is(PD,:);
   PDremap = zeros(size(subremap),useInt);
   PDremap(keep) = ismi+ismiconj;
   %% Region that is the complex conjugate of the principal domain
   PDconjugate = false(size(subremap));
%    PDconjugate(unaliased) = sum(sign(Wsrt(1:3,:)))<0;
   PDconjugate(keep) = ismconj;
   
%   Isreduced = IsreducedPD;
   Is = IsPD;
   subremap = PDremap;
   
   %%% Identify the symmetry type for each region. At least one region
   %%% from each group must be estimated for the partial delay estimator.
   %%% Symmtery types correspond to permutations of position of the first
   %%% term in the average cross-polyspectral product: a point is in region
   %%% k when its frequencies 2..order are, as a set, the frequencies of its
   %%% principal-domain representative without the k-th one (or their
   %%% negation), i.e. when the point's first frequency sits in the k-th
   %%% slot of the representative; the smallest such k is taken. Matching
   %%% against the point's own representative is the same as matching
   %%% against every principal-domain tuple, since a match fixes the tuple
   %%% up to sign.
   W23 = W(2:order,:);
   [W23srt,w23srti] = sort(round(abs(W23)./tol)*tol);
   w23srti = w23srti + (order-1)*repmat(0:size(W23,2)-1,order-1,1);
   W23srt = W23srt.*sign(W23(w23srti));
    SR = uint8(zeros(sum(keep(:)),1));
    rep = ismi+ismiconj;                        % representative of each kept point (0 = none)
    has = rep > 0;
    WPD = W(:,PD);
    canon = @(X) canonicalForm(X,tol);
    W23has = W23srt(:,has);
    for k = order:-1:1
        W23pdall = canon(WPD(setdiff(1:order,k),:));      % canonical minus-k form of every PD tuple
        W23pdsrt = W23pdall(:,rep(has));                  % ... of each point's representative
        sri = all(W23has == W23pdsrt,1) | all(W23has == -W23pdsrt,1);
        srk = SR(has);
        srk(sri) = k;
        SR(has) = srk;
        if any(~has)
            %%% A point without a representative (the first kept point is
            %%% never a conjugate, see above) is matched against every
            %%% principal-domain tuple, as all points were before.
            sri0 = ismember(W23srt(:,~has)',[W23pdall';-W23pdall'],'rows');
            srk = SR(~has);
            srk(sri0) = k;
            SR(~has) = srk;
        end
    end
 
    
 

else 
%    PDconjugate = zeros(size(remap));
%    PDconjugate(unaliased) = sum(sign(Wsrt(1:3,:)))<0;
   
    PDconjugate = false;
end



% keeplp = arrayfun(@(fr,lp)abs(fr{1})<=lp,freqsin(1:end-1),lowpass(1:end-1),'uniformoutput',false);
% keeplp2 = cellfun(@(kpfr,kplp)kpfr(kplp),keepfreqs(1:end-1),keeplp,'uniformoutput',false);
% dims = cellfun(@(x)sum(x),keeplp);

% keepregion = false(dims);
% keepregion(keeplp2{:})=true;
keepall = false(dims);
keepall(keepregion) = keep;
Z = zeros(dims,useInt);
PDall=Z;
PDall(keepall)=PD;
remap = Z;
remap(keepregion) = subremap;

remap(remap==0) = size(Is,1)+1;
out.Is = Is;%Isreduced;
out.freqs = freqs(1:end-1);
out.keep = keepall;
out.principal_domain = PDall;
out.remap = remap;
out.reduce = find(PDall);
PDconj = false(size(Z));
PDconj(keepregion) = PDconjugate;
out.PDconj = PDconj;
out.Bfreqs = cellfun(@(fr,kpfr)fr(kpfr),freqsin(1:end-1),keeplp,'uniformoutput',false);
SymReg = uint8(zeros(size(keepall)));
SymReg(keepall)=SR;
out.partialSymmetryRegions=SymReg;

function [tf,loc] = rowsMember(A,B,tol)

% [tf,loc] = rowsMember(A,B,tol): ismember(A',B','rows') for canonical
% tuples whose entries are signed multiples of tol, through one integer
% key per column when the keys fit in a uint64 (exact), otherwise through
% ismember itself.

order = size(A,1);
Q = round(abs([A B])./tol);
sgn = [A B] < 0;
base = 2*max(Q(:))+2;                          % codes 2*Q + sign lie in [0, base-1]
if isempty(Q) || base^order < 2^63
    code = uint64(2*Q + sgn);
    key = zeros(1,size(code,2),'uint64');
    b = uint64(1);
    for k = 1:order
        key = key + code(k,:).*b;
        b = b.*uint64(base);
    end
    [tf,loc] = ismember(key(1:size(A,2)),key(size(A,2)+1:end));
    tf = tf(:); loc = loc(:);
else
    [tf,loc] = ismember(A',B','rows');
end

function C = canonicalForm(X,tol)

% Columns of X rounded to tol, sorted by magnitude, with their signs (the
% form the symmetry-region comparison uses; no tie-break between equal
% magnitudes of opposite sign, as before).

[C,srti] = sort(round(abs(X)./tol)*tol);
srti = srti + size(X,1)*repmat(0:size(X,2)-1,size(X,1),1);
C = C.*sign(X(srti));
