
function A = pica(X0,ncomp,ord,a0,dnl)

%ICA through power iteration

if nargin < 2 || isempty(ncomp)
    ncomp = size(X0,2);
end

if nargin < 3 || isempty(ord)
    ord = 4;
end

if nargin < 4 || isempty(a0)
    a0 = randn(size(X0,2),1);
end

maxiter = 500;
%%

m0 = mean(X0);
R0 = cov(X0);
WM = R0^(.5);
Z0 = (X0-m0)*pinv(WM);

X = Z0;

if isnumeric(ord)
    nl = @(x)x.^ord;
    dnl = @(x)ord*x.^(ord-1);
elseif ischar(ord)
    switch ord
        case 'logcosh'
            nl = @(x)log(cosh(x));
            dnl = @(x)tanh(x);
    end
end

for dim = 1:ncomp

R = X'*X;
Rinv = pinv(R);



tol = 1e-6;
d = Inf;

a = WM*a0;

% cm = cumulant(X*a,ord);

iter = 1;
clear ds kt sk
while d> tol && iter < maxiter
    
    r = X*a;


    if isnumeric(ord)
         rX = r.^((ord-2)/2).*X; 
         G = (rX'*rX)*Rinv;
    else
        rX = sign(r).*dnl(r).*X./(abs(r) +eps) ;
        G = (rX'*X)*Rinv;
    end
    anew = G'*a;
%     
     anew = anew./norm(anew);
    d = norm(anew-a);
    
%     cm(iter+1) = cumulant(X*anew,4);
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
A = A.*sign(mean(nl(X0*A)));
