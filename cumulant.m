function cum = cumulant(X,order,dim,normalize)

% cum = cumulant(X,order,dim,normalize)
%Compute the order-th sample cumulant of X along dimension dim (default 1),
%ignoring NaN. With normalize (default true) the cumulant is divided by the
%variance to the power order/2: the standardized cumulant (the skewness for
%order 3, the excess kurtosis for order 4).
%
%Each raw moment mu_j = nanmean(X.^j,dim), j = 1..max(order,2), is one pass
%over X; the cumulants follow from the moment recursion
%   kappa_m = mu_m - sum_{k=1}^{m-1} nchoosek(m-1,k-1)*kappa_k*mu_{m-k}.

if nargin < 3 || isempty(dim)
    dim = 1;
end

if nargin < 4 || isempty(normalize)
    normalize = true;
end

% The operations, and their order, are those of the self-recursive r1811
% version, so for a double order along dim 1 of an array with more than one
% row the results are bit-identical to it: keep X.^j (not repeated .*X),
% nanmean, and k running upward.
n = max(order,2);
mu = cell(1,n);
for j = 1:n
    mu{j} = nanmean(X.^j,dim);
end

kap = mu;
for m = 2:n
    for k = 1:m-1
        kap{m} = kap{m} - nchoosek(m-1,k-1)*kap{k}.*mu{m-k};
    end
end

cum = kap{order};
if normalize
    cum = cum./kap{2}.^(order/2);
end
