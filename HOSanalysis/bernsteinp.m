function p = bernsteinp(n,order)

if isscalar(n)
    x = linspace(0,1,n);
else
    x = n;
end
for k = 0:order
    
    p(:,k+1) = nchoosek(order,k)*x.^k.*(1-x).^(order-k);
end