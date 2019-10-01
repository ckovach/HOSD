function [permP,Ntot] = phase_permutation_test(hos,x,maxperm)

%
% [permP,ntot]= phase_permutation_test(hos, x, maxperm)
%
% Permutation test on bispectral phase. 
% At each permutation, the phase of the deterministic bispectrum for each
% window in the estimator is randomized and the the result is averaged to
% obtain a surrogate sample bispecral estimate.
%
% For the sake of efficiency, the kth permutation is computed only at
% coefficients for which the permutation P-value lies within 2 standard
% errors of the Bonferroni correction threshold and continues until all
% coefficients meet the stopping criterion or until maxperm is reached.
%
%    Inputs:
%
%  hos : hosobject used to compute the bispectrum
%  x   : data in the form of a column vector or N x M matrix
%         where N = hos.buffersize
%  maxperm : maximum number of permutations (default = 25e3)
%          
%   Outputs: 
%
%  permP : Permutation p-value. This has a minimum value of 0.5/maxperm.
%  ntot  : Total number of permutations conducted for each coefficient.
%
% See HOSOBJECT

%C. Kovach 2019
alpha = .05;
if nargin < 3 || isempty(maxperm)
    maxperm = 25e3; %%% Maximum number of permutations. If set to Inf, then will
                   %%% continue until all P-values are 2 std errors outsid
                   %%% of the bonferonni threshold.
end

stop_threshold = 2; %%% For the sake of efficiency, exclude coefficients after
                        %%% permutation P value differs from the bonferroni threshold by this many
                        %%% std. errors.

if size(x,1)~=hos(1).buffersize
    X = hos(1).chop_input(x);
end
isn = any(isnan(X));

if sum(isn)>0
    fprintf('\n%i windows (%0.2f%%) containing nans will be ignored.',sum(isn),100*sum(isn)/length(isn))
    X(:,isn)=[];
end

FX = fft(X);

%%% For the purpose of this test, we only care if the numerator differs
%%% significantly from zero, hence no need for normalization. This is
%%% because the denominator does not depend on phase and so does not change
%%% across permutations. The result applies both to the bispectral estimtae
%%% and bicoherence.

FFX = conj(FX(hos(1).freqindx.Is(:,hos(1).order),:));
        
for k = hos(1).order-1:-1:1  
   FFX = FFX.*FX(hos(1).freqindx.Is(:,k),:);
end

B0 = abs(mean(FFX,2));

pbonf =alpha./size(B0,1); %Bonferroni threshold will be used in the stopping criterion.
                        %Permutations will continue only for points that
                        %are within 3 std err of the Bonferroni threshold
                        %(until maxperm is reached).
                        
nsig = zeros(size(B0,1),1);
ntot = nsig;
keep = true(size(nsig));

permi = 1;
fpn = 0;

reseed

keep0=true;
while any(keep) &&  permi<maxperm

    Bperm = mean(abs(FFX(keep,:)).*exp(2*pi*1i.*rand(size(FFX(keep,:)))),2);
    
    nsig(keep) = nsig(keep)+(abs(Bperm)>=B0(keep));
    
    ntot = ntot+keep;
    
    pperm = (nsig+.5)./(ntot+1); % Regularize the pvalue estimate so that it is never 1 or 0.
    
    keep = abs((pperm-pbonf).*sqrt(ntot./(pperm.*(1-pperm))))<stop_threshold;
%    keep = nsig<stop_threshold;
    fpn = fprintf([repmat('\b',1,fpn),'\nperm %i (N remaining = %i)'],permi,sum(keep))-fpn;
    
    permi = permi+1;
    
    keep0=keep;
end
    
permP = pperm;
permP(end+1) = nan;
permP = permP(hos(1).fullmap);

if nargout > 1
    Ntot = ntot;
    Ntot(end+1) = 0;
    Ntot = Ntot(hos(1).fullmap);
end