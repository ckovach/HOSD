function cum = cumcum(X,order,dim)

% cum = cumulant(X,n)
%Compute a cumulative nth order cumulant along the first dimension of X.

if nargin < 3 || isempty(dim)
    dim = 1;
end

cmlN = cumsum(ones(size(X)),dim);

cum = cumsum(X.^order)./cmlN;

for k = 1:order-1
   
    cum = cum - nchoosek(order-1,k-1)*cumcum(X,k,dim).*(cumsum(X.^(order-k))./cmlN);
    
end