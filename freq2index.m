function out = freq2index(freqs,order,lowpass,highpass,freqindex,condense,frequency_spacing)

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
    if iscell(freqs)
        order = length(freqs);
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

%%
if isnumeric(freqs)
    freqs = repmat({freqs},1,order);
    if nargin < 6 || isempty(condense)
        condense = true;
    end
elseif nargin < 6 || isempty(condense)
    condense = false; % If the indexing is the same for the different frequencies, then we will give the mapping into the principal domain
end
if nargin < 5 || isempty(freqindex)
   freqindex = cellfun(@(x)1:length(x),freqs,'uniformoutput',false);
elseif isnumeric(freqindex)
    freqindex = repmat({freqindex},1,order);
elseif length(freqindex)<order
    freqindex(end+1:order) = freqindex(end);
end

%%
%%% For cross-polyspectra involving fewer signals than the specified order,
%%% assume that the last signal is repeated.
if length(freqs) < order
    freqs(end+1:order) = freqs(end);
end

%%% Check if negative frequencies are explicitly represented
%%% Will index into the negative frequencies if so.
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

frsrti(end+1) = length(freqs{order})+1;
[srt,srti] = sort([(-1)^two_sided*Fsum(:);frcent(:)]);
E = [zeros(n,1);ones(length(frcent),1)];
E = E(srti); 
IND = zeros(size(E));
IND(srti) = cumsum(E);
IND(n+1:end) = [];
IND(IND==0) = length(freqs{end})+1;
IND =  frsrti(reshape(IND,size(Fsum)));

Is{order}=IND(:) ;

% keep = IND>0 & IND<=length(freqs{end});


Is = cellfun(@(x,fri)fri(x(:))',Is,freqindex,'uniformoutput',false);
Is = [Is{:}];
%Isreduced = Is(keep,:);
W = [];
for k = 1:order
   W(k,:) = Ws{k};%(keep);
end
remap = zeros(size(keep));
remap(keep) = find(keep);

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
   PDremap = zeros(size(remap));
   PDremap(keep) = ismi+ismiconj;
   %% Region that is the complex conjugate of the principal domain
   PDconjugate = false(size(remap));
%    PDconjugate(unaliased) = sum(sign(Wsrt(1:3,:)))<0;
   PDconjugate(keep) = ismconj;
   
%   Isreduced = IsreducedPD;
   Is = IsPD;
   remap = PDremap;
else 
%    PDconjugate = zeros(size(remap));
%    PDconjugate(unaliased) = sum(sign(Wsrt(1:3,:)))<0;
   
    PDconjugate = false;
end

remap(remap==0) = size(Is,1)+1;
out.Is = Is;%Isreduced;
out.freqs = freqs(1:end-1);
out.keep = keep;
out.principal_domain = PD;
out.remap = remap;
out.PDconj = PDconjugate;