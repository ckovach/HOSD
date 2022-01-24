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

% Copyright Christopher Kovach, University of Iowa, 2018

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



if  diagonal_slice
    
    dims = cellfun(@length,freqs(1:order-1));
    
    dgs = sparse(prod(dims),1);
    repfr =[1 cumprod(dims)];
    keep = false(dims);
    ws = repmat({[]},1,order-1);
    [ws{:}] = ndgrid(freqs{1:order-1});
    word = -sum(cat(order,ws{:}),order);
%     dgs(:) = sum(abs(cat(order,ws{:})-word)<eps.*abs(word),order);
    dgtol = 1.5*me.sampling_rate/me.buffersize;
    dgs(:) = sum(abs(cat(order,ws{:})-word)< dgtol,order);
 
    for k = 1:order-1
     
        for kk = k+1:length(freqs)-1
           
            iseq = sparse(freqs{k}'==freqs{kk});
            
            rpm=repmat(kron(iseq,ones(repfr(k),1)),repfr(kk-1)./repfr(k),1);
            rec = repmat(rpm(:),prod(dims)/numel(rpm),1);
%             rec = kron(ones(prod(dims)/numel(rpm)),rpm(:));
            
            dgs = dgs+rec(:);
%             dgs(:) = dgs(:) + (ws{k}(:)==ws{kk}(:));
        end
%         dgs = kron(ones(numel(iseq),1),dgs)+iseq(:);
    end
    
    keep(:) =dgs>=diagonal_slice;   
    
    Is = repmat({[]},1,order-1);
    [Is{:}] = ind2sub(dims,find(keep));
    kpindx = keep(:);
else
%     keep  = true;
    
    freqsi = cellfun(@(x)1:length(x),freqs(1:end-1),'uniformoutput',false);
    Is = repmat({[]},1,order-1);
    [Is{:}] = ndgrid(freqsi{:});
    Is = cellfun(@(x)x(:),Is,'uniformoutput',false);
    keep = true(size(Is{1}(:)));
    kpindx = keep;
end


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
    Wseq = 0;
    for k = 1:length(Ws)-1
        for kk = k+1:length(Ws)
           Wseq = Wseq + (abs(Ws{k} - Ws{kk}) < dgtol);
        end
    end

    keep(keep) = Wseq>=diagonal_slice;      
% else
%     keep  = true;
end

if any(~isinf(xlowpass))
    keepxlp=false;
    for k = 1:length(Ws)
        keepxlp = keepxlp | abs(Ws{k})<xlowpass(k);
    end
     keep(kpindx) = keepxlp & keep(kpindx);
%     keep = keepxlp & keep;
end

if any(xhighpass>0)
    keepxhp=false;
    for k = 1:length(Ws)
        keepxhp = keepxhp | abs(Ws{k})>=xhighpass(k);
    end
     keep(kpindx)  = keep(kpindx) & keepxhp;
%     keep  = keep & keepxhp;
end

passtol = .25*me.sampling_rate/me.buffersize;
for k = 1:length(Ws)
   
    if ~isinf(lowpass(k))
         keep(kpindx) = keep(kpindx) & abs(Ws{k})<lowpass(k)-passtol;
%         keep = keep & abs(Ws{k})<lowpass(k);
    end
    if highpass(k)>0
%         keep = keep & abs(Ws{k})>highpass(k);
         keep(kpindx) = keep(kpindx) & abs(Ws{k})>highpass(k)+passtol;
    end
end
if ~isscalar(mask)
    keep(:) = keep(:)&mask(:);
end

if ~isempty(include_sig)
    Wsig=abs(sum(sign([Ws{:}]),2));
    
    keep(kpindx) = keep(kpindx) & ismember(Wsig,include_sig);    
end

for k = 1:length(Ws)
     Ws{k} = Ws{k}(keep(kpindx));
%     Ws{k} = Ws{k}(keep);
    if k<=length(Is)
         Is{k} = Is{k}(keep(kpindx));
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
        discard = false(size(keep));
        WKP = [Ws{:}];
     %   WKP = WKP(keep,:);
        for k = 1:order-1
            discard(keep) = discard(keep) | any(abs(repmat(WKP(:,k),1,order-k)+WKP(:,k+1:order)) < max(shighpass)+passtol,2)   | all(abs(repmat(WKP(:,k),1,order)+WKP) > max(slowpass)-passtol,2)  ;       
        end

        PD(discard(keep)) = [];
        for k = 1:length(Ws)
            Ws{k} = Ws{k}(~discard(keep));
            if k<=length(Is)
                Is{k} = Is{k}(~discard(keep));
            end
        end
        keep(keep) = keep(keep) & ~discard(keep);
    end
end





 
