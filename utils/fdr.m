function [fdrout, q, pi0, rs] = mafdr(p, varargin)

% Compute the fdr using BH Linear step-up procedure
% BHFDR taken from the matlab function MAFDR

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

nanp = isnan(p);

p = p(~nanp);
m = numel(p);
[p_ord, idx] = sort(p);

% Find cutoffs starts from end
fdr_ord = zeros(m,1, class(p));
fdr_ord(m) = p_ord(m);

for k = m-1:-1:1
    fdr_ord(k) = p_ord(k)*m/k;
    if fdr_ord(k) > fdr_ord(k+1)
        fdr_ord(k) = fdr_ord(k+1);
    end
end
fdr(idx) = fdr_ord;
fdrout(~nanp) = fdr(:);
fdrout(nanp) = nan;




    


