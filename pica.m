
function A = pica(X0,ord,a0)

%ICA through power iteration

if nargin < 2 || isempty(ord)
    ord = 4;
end

if nargin < 3 || isempty(A0)
    a0 = randn(size(X0,2),1);
end

maxiter = 500;
%%

m0 = mean(X0);
R0 = cov(X0);
WM = R0^(.5);
Z0 = (X0-m0)*pinv(WM);

X = Z0;



for dim = 1:size(X,2)

R = X'*X;
Rinv = pinv(R);



tol = 1e-6;
d = Inf;

a = a0;

cm = cumulant(X*a,ord);

iter = 1;
clear ds kt sk
while d> tol && iter < maxiter
    
    r = X*a;

%     rX = r.^(ord-2).*X; 
%     G = (rX'*X)*Rinv;
     rX = r.^((ord-2)/2).*X; 
     G = (rX'*rX)*Rinv;
    
%       [anew,ev(iter)] = eigs(G,1);
    anew = G'*a;
%     
     anew = anew./norm(anew);
    d = norm(anew-a);
    
%     cm(iter+1) = cumulant(X*anew,ord);
%     if cm(iter+1)>cm(iter)
        a = anew;
%     else
%         d=-Inf;
%     end
%     ds(iter)=d;
%     kt(iter) = kurtosis(r);
%     sk(iter) = skewness(r);
   iter = iter+1;
end
A(:,dim) = a;   
%%
 X = X*(eye(size(X,2))-a*a');
end
A = pinv(WM)*A;
