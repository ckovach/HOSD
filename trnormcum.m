function c = trnormcum(a,b,ord)

%Cumulants of the truncated normal distribution

if isinf(a)
    a = sign(a)*100;
end
if isinf(b)
    b = sign(b)*100;
end

m(1,length(b)) = 0;
m(2,1:length(b)) = 1;
for k = 1:ord

    %Orjeban 2014
    m(k+2,:) = (a.^(k-1).*normpdf(a) - b.^(k-1).*normpdf(b))./(normcdf(b)-normcdf(a)) + (k-1)*m(k,:);
    
end

    c = m(end,:);

for k = 1:ord-1
   c = c -  nchoosek(ord-1,k-1)*trnormcum(a,b,k).*m(ord-k+2,:);
end
