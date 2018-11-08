function out = freq2index(freqsin,order,lowpass,highpass,keepfreqs,condense,frequency_spacing)

% [Is,remap] = freq2index(freqs,order)
%
% Function to generate indexing for polyspectra of a given order.
%
% INPUT:
%
% freqs - vector or 1 x order cell array of vectors containing frequencies.
%
%

% Copyright Christohpher Kovach, University of Iowa 2018.

if nargin < 7 || isempty(frequency_spacing)
    frequency_spacing = 'linear';
end
if nargin < 2 || isempty(order)
    if iscell(freqsin)
        order = length(freqsin);
    else 
        order = 3;
    end
end
if nargin < 3 || isempty(lowpass)
    lowpass = .5;
end
if nargin < 4 || isempty(highpass)
    highpass = 0;
end

if isscalar(lowpass)
    lowpass = ones(1,order)*lowpass;
end

if isscalar(highpass)
    highpass= ones(1,order)*highpass;
end
%%
if isnumeric(freqsin)
    freqsin = repmat({freqsin},1,order);
    if nargin < 6 || isempty(condense)
        condense = true;
    end
elseif nargin < 6 || isempty(condense)
    condense = false; % If the indexing is the same for the different frequencies, then we will give the mapping into the principal domain
end
if nargin < 5 || isempty(keepfreqs)
   keepfreqs = cellfun(@(x)1:length(x),freqsin,'uniformoutput',false);
elseif isnumeric(keepfreqs)
    keepfreqs = repmat({keepfreqs},1,order);
elseif length(keepfreqs)<order
    keepfreqs(end+1:order) = keepfreqs(end);
end

%%
%%% For cross-polyspectra involving fewer signals than the specified order,
%%% assume that the last signal is repeated.
if length(freqsin) < order
    freqsin(end+1:order) = freqsin(end);
end

%%% Check if negative frequencies are explicitly represented
%%% Will index into the negative frequencies if so.
freqs = cellfun(@(fr,kp)fr(kp),freqsin,keepfreqs,'uniformoutput',false);

nneg = sum(freqs{end}<0);
npos = sum(freqs{end}>0);
two_sided = nneg>npos/2; %%% Assume that negative frequencies are just padding if there aren't at least as many negative frequencies as half the number of positive frequencies.

[PD,Ws,Is,keep] = find_principal_domain(freqs,order,lowpass,highpass);
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
[srt,srti] = sort([(-1)^two_sided*Fsum(:);frcent(:)]);
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

subremap = zeros(size(keep));
subremap(keep) = find(keep);

tol = min(abs(diff([freqs{:}])))/2;

if condense %for auto-spectra we only need the principal domain. This is not so for cross-spectra
    [Wsrt,wsrti] = sort(round(abs(W)./tol)*tol);
    wsrti = wsrti + order*repmat(0:size(W,2)-1,order,1);
    Wsrt = Wsrt.*sign(W(wsrti));
   
     [~,ismi] = ismember(Wsrt',Wsrt(:,PD)','rows');
   [ismconj,ismiconj] = ismember(-Wsrt',Wsrt(:,PD)','rows');
   ismiconj(1) = 0;
   ismconj(1)=false;
%   IsreducedPD = Isreduced(PD(keep),:);
   IsPD = Is(PD,:);
   PDremap = zeros(size(subremap));
   PDremap(keep) = ismi+ismiconj;
   %% Region that is the complex conjugate of the principal domain
   PDconjugate = false(size(subremap));
%    PDconjugate(unaliased) = sum(sign(Wsrt(1:3,:)))<0;
   PDconjugate(keep) = ismconj;
   
%   Isreduced = IsreducedPD;
   Is = IsPD;
   subremap = PDremap;
   
   %%% Identify the symmetry regions for partial cross bicoherence
   W23 = W(2:order,:);
      [W23srt,w23srti] = sort(round(abs(W23)./tol)*tol);
    w23srti = w23srti + (order-1)*repmat(0:size(W23,2)-1,order-1,1);
    W23srt = W23srt.*sign(W23(w23srti));
    SR = uint8(zeros(sum(keep(:)),1));
    for k = order:-1:1
        W23pd=W(setdiff(1:order,k),PD);
        [W23pdsrt,w23pdsrti] = sort(round(abs(W23pd)./tol)*tol);
        w23pdsrti = w23pdsrti + (order-1)*repmat(0:size(W23pdsrt,2)-1,order-1,1);
        W23pdsrt = W23pdsrt.*sign(W23pd(w23pdsrti));
        sri = ismember(W23srt',[W23pdsrt';-W23pdsrt'],'rows');
        SR(sri) = k;
    end
 
    
 

else 
%    PDconjugate = zeros(size(remap));
%    PDconjugate(unaliased) = sum(sign(Wsrt(1:3,:)))<0;
   
    PDconjugate = false;
end

keeplp = arrayfun(@(fr,lp)abs(fr{1})<=lp,freqsin(1:end-1),lowpass(1:end-1),'uniformoutput',false);
keeplp2 = cellfun(@(kpfr,kplp)kpfr(kplp),keepfreqs(1:end-1),keeplp,'uniformoutput',false);
dims = cellfun(@(x)sum(x),keeplp);

keepregion = false(dims);
keepregion(keeplp2{:})=true;
keepall = false(dims);
keepall(keepregion) = keep;
Z = zeros(dims);
PDall=Z;
PDall(keepall)=PD;
remap = Z;
remap(keepregion) = subremap;

remap(remap==0) = size(Is,1)+1;
out.Is = Is;%Isreduced;
out.freqs = freqs(1:end-1);
out.keep = keep;
out.principal_domain = PD;
out.remap = remap;
out.reduce = find(PDall);
PDconj = false(size(Z));
PDconj(keepregion) = PDconjugate;
out.PDconj = PDconj;
out.Bfreqs = cellfun(@(fr,kpfr)fr(kpfr),freqsin(1:end-1),keeplp,'uniformoutput',false);
SymReg = uint8(zeros(size(keepall)));
SymReg(keepall)=SR;
out.partialSymmetryRegions=SymReg;
