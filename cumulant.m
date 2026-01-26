function cum = cumulant(X,order,dim,normalize)

% cum = cumulant(X,n)
%Compute the nth order sample cumulant along the first dimension of X.

if nargin < 3 || isempty(dim)
    dim = 1;
end

if nargin < 4 || isempty(normalize)
    normalize = true;
end

cum = nanmean(X.^order,dim);

for k = 1:order-1
   
    cum = cum - nchoosek(order-1,k-1)*cumulant(X,k,dim,false).*nanmean(X.^(order-k));
    
end

if normalize
    cum = cum./cumulant(X,2,dim,false).^(order/2);
end
