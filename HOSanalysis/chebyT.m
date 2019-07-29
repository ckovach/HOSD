function P = chebyT(n,order)

% P = chebyT(n,order)
%Return Chebyshev polynomials up to the given order for a window of size n.


P = zeros(n,order+1);

P(:,1)=1;
if order > 0
    P(:,2) = linspace(-1,1,n);
else
    return
end

for k = 3:order+1
   
    P(:,k) = 2*P(:,2).*P(:,k-1) - P(:,k-2);
    
end