function [permP,Ntot,Q] = phase_permutation_test(hos,x,maxperm,stop_threshold)

%
% [permP,ntot]= phase_permutation_test(hos, x, [maxperm],[stop_threshold)
%
% Permutation test on bispectral phase. 
% At each permutation, the phase of the deterministic bispectrum for each
% window in the estimator is randomized and the the result is averaged to
% obtain a surrogate sample bispectral estimate.
%
% For the sake of efficiency, the kth permutation is computed only at
% coefficients for which the current permutation P-value lies within 2 standard
% errors of the Bonferroni correction threshold and continues until all
% coefficients meet the stopping criterion or until maxperm is reached.
%
%    Inputs:
%
%  hos : hosobject used to compute the bispectrum
%  x   : data in the form of a column vector or N x M matrix
%         where N = hos.buffersize
%  maxperm : maximum number of permutations (default = 25e3)
%  stop_threshold:
%            if 0<  . <1     :        the p-value threshold used in place of
%                                     bonferroni for stopping.
%            if . > 1        :        the number of std errors from the Bonferroni
%                                     threshold to tolerate. 
%            if . =[pthresh, sdetol] : first value is p threshold, second is
%                                     std err tolerance. 
%   Outputs: 
%
%  permP : Permutation p-value. This has a minimum value of 0.5/maxperm.
%  ntot  : Total number of permutations conducted for each coefficient.
%
% See HOSOBJECT

%C. Kovach 2019
alpha = .05;
fdrthresh = .01;
default_stop_threshold = 2;
use_fdr_thresh = true; %%% Adjusts stopping threshold according to false discovery rather than Bonferroni correction, if true.
                       %%% This will often be much more efficient.
                       
if nargin < 3 || isempty(maxperm)
    maxperm = 1e6; %%% Maximum number of permutations. If set to Inf, then will
                   %%% continue until all P-values are 2 std errors away
                   %%% from the Bonferonni or FDR threshold.
end

if nargin < 4 || isempty(stop_threshold)
    stop_threshold = default_stop_threshold; %%% For the sake of efficiency, exclude coefficients after
                        %%% permutation P value differs from the Bonferroni or FDR threshold by this many
                        %%% std. errors.
end
    
if size(x,1)~=hos(1).buffersize
    X = hos(1).chop_input(x);
else
    X = x;
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

B0 = abs(mean(FFX,2)); %%% Magnitude of the bispectral estimate

if stop_threshold(1) < 1
    if use_fdr_thresh
        fdrthresh = stop_threshold(1);
    else
        pbonf = stop_threshold;
    end
    if length(stop_threshold)>1
        stop_threshold(2)=1;
    else
        stop_threshold=default_stop_threshold;
    end 
else
    pbonf =alpha./size(B0,1); %Bonferroni threshold will be used in the stopping criterion.
end                        
                        
nsig = zeros(size(B0,1),1);
ntot = nsig;
keep = true(size(nsig));

permi = 1;
fpn = 0;

reseed

while any(keep) &&  permi<maxperm 
    
    %%% Surrogate estimate for which phase has been randomized
    Bperm = mean(abs(FFX(keep,:)).*exp(2*pi*1i.*rand(size(FFX(keep,:)))),2);
    
    %%% Number of surrogate sample estimates with magnitude equal or greater than the original estimate
    nsig(keep) = nsig(keep)+(abs(Bperm)>=B0(keep));
    
    %%% Total number of permutations
    ntot = ntot+keep;
    
    %%% Regularized permutation p-value
    pperm = (nsig+.5)./(ntot+1); % Regularize the pvalue estimate so that it is never 1 or 0.
                                 % This is done here by adding 0.5 to the
                                 % numerator and 1 to the denominator.
                                 % The minimum p-value will therefore be
                                 % 0.5/(maxperm + 1).
                                 
    %%% Find coefficients that require no further testing.
    keep = abs((pperm-pbonf).*sqrt(ntot./(pperm.*(1-pperm))))<stop_threshold;
%    keep = nsig<stop_threshold;
    
    if mod(permi,ceil(.001*permi)*10) ==0
        fpn = fprintf([repmat('\b',1,fpn),'\nperm %i (N remaining = %i, pthresh = %0.2g)'],permi,sum(keep),pbonf)-fpn;
        if use_fdr_thresh 
            %%%Adjust the threshold according to the false discovery
            %%%criterion (Benjamini-Hochberg)
           srtp = sort(pperm); 
           srtp(1)=0;
           pbonf = find(srtp'<(1:length(srtp))/length(srtp)*fdrthresh,1,'last')/length(pperm)*fdrthresh;
        end
        if sum(keep)/(sum(pperm<pbonf & ~keep)+1)<.05 % If the number of tests within error range of the stopping threshold
                                                      % is less than 5% of tests deemed significant, then stop. 
            maxperm=0;
        elseif (sqrt(pbonf*(1-pbonf))./sqrt(permi))/pbonf < .1 %If error range is less than 10% of the threshold, then stop.
            maxperm = 0;
        end
    end
    permi = permi+1;

end

%%% Reshape 
permP = pperm;
permP(end+1) = nan;
permP = permP(hos(1).fullmap);

if nargout > 1
    Ntot = ntot;
    Ntot(end+1) = 0;
    Ntot = Ntot(hos(1).fullmap);
end
if nargout > 2
    q = fdr(pperm);
    q(end+1)=nan;
    Q = q(hos(1).fullmap);
end