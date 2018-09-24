function [PD,signatures,Ws,Is] = find_principal_domain(freqs,order)

% Find the principal domain in a higher-order spectrum
%
% The principal domain represents a set of non-redundant regions in the
% higher order spectrum. Redundancy has two sources: (1)
% permutation of axes and (2) conjugate symmetry. By (1) we can include
% only regions such that W1<W2<...<Wk<0 and by (2) we 
% can exclude regions with signature given by a sign reversal of included region. 
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

if nargin < 2 || isempty(order)
    if iscell(freqs)
        order = length(freqs);
    else 
        order = 3;
    end
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


%%% The number of regions in the principal domain depends on the possible 
%%% signatures.
nsig = floor(order/2);
signatures = (-1).^(repmat(1:order,nsig,1)- 1>= order-repmat((1:nsig)',1,order));

freqsi = cellfun(@(x)1:length(x),freqs(1:end-1),'uniformoutput',false);
Is = repmat({[]},1,order-1);
[Is{:}] = ndgrid(freqsi{:});

Ws = {};
Wsum=0;
for k = 1:length(Is)
    Ws{k} = freqs{k}(Is{k});
    Wsum  = Wsum+Ws{k};
end
Ws{order} = -Wsum;
    
PD = false;
for k = 1:nsig
      PD0 = Ws{1}>=0 & Ws{order}<=0; %First signature is always + and last always -.    
         for kk = 2:order-1
             PD0 = PD0 & signatures(k,kk)*Ws{kk}>=signatures(k,kk-1)*Ws{kk-1};
         end
 

    PD = PD | PD0;
  
end
 