function out = freq2index(freqs,order,frequency_spacing)

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

if nargin < 3 || isempty(frequency_spacing)
    frequency_spacing = 'linear';
end
if nargin < 2 || isempty(order)
    if iscell(freqs)
        order = length(freqs);
    else 
        order = 3;
    end
end
    
%%
condense = false; % If the indexing is the same for the different frequencies, then we will give the mapping into the principal domain
if isnumeric(freqs)    
    freqs = repmat({freqs},1,order);
    condense = true;
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

[PD,signatures,Ws,Is] = find_principal_domain(freqs);
Fsum = Ws{end};
% frinds = cellfun(@(x)1:length(x),freqs,'uniformoutput',false);
% 
% Is=repmat({[]},1,order-1);
% [Is{:}] = ndgrid(frinds{1:order-1});
% 
% Fsum = 0;
% for k =1:order-1
%     Fsum = Fsum + freqs{k}(Is{k});
% end

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

Is{order}=IND ;

unaliased = IND>0 & IND<=length(freqs{end});


Is = cellfun(@(x)x(:),Is,'uniformoutput',false);
Is = [Is{:}];
Isreduced = Is(unaliased,:);
W = [];
for k = 1:order
   W(k,:) = Ws{k}(unaliased);
end
remap = zeros(size(IND));
remap(unaliased) = find(unaliased);
% CPconj(unaliased,order-1) = false;
% CPconj(unaliased,order) = ~two_sided;
%
% Determine principal domain
% PD = unaliased;
% PD(unaliased)= sum(W(1:order-1,:)) + 2*W(order,:)>0;
% for k = 2:size(W,1)
%     PD(unaliased) = PD(unaliased) & W(order,k-1)>=W(order,k);
% end
% 

tol = min(abs(diff([freqs{:}])))/2;
if condense
   
   [Wsrt,wsrti] = sort(round(abs(W)./tol)*tol);
   wsrti = wsrti + order*repmat(0:size(W,2)-1,order,1);
   Wsrt = Wsrt.*sign(W(wsrti));
   [ism,ismi] = ismember(Wsrt',Wsrt(:,PD(unaliased))','rows');
   [ismconj,ismiconj] = ismember(-Wsrt',Wsrt(:,PD(unaliased))','rows');
   ismiconj(1) = 0;
   ismconj(1)=0;
   IsreducedPD = Isreduced(PD(unaliased),:);
   PDremap = zeros(size(remap));
   PDremap(unaliased) = ismi+ismiconj;
   %% Region that is the complex conjugate of the principal domain
   PDconjugate = zeros(size(remap));
   PDconjugate(unaliased) = ismconj;
   Isreduced = IsreducedPD;
   remap = PDremap;
else 
    PDconjugate = false;
end

remap(remap==0) = size(Isreduced,1)+1;
out.Is = Isreduced;
out.freqs = freqs(1:end-1);
out.unaliased = unaliased;
out.principal_domain = PD;
out.remap = remap;
out.PDconj = PDconjugate;
