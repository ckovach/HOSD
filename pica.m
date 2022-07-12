
function [UM,WM,A] = pica(X0,ncomp,ord,a0,dnl,verbose)

%[UM,WM,A] = pica(X0,ncomp,ord,a0,dnl,verbose)
%ICA through power iteration. This finds a series of projection that
%maximize the moment of a given order. By default the 4th moment is used, which 
%assumes all independent components are leptokurtic.

if nargin < 2 || isempty(ncomp)
    ncomp = size(X0,2);
end

if nargin < 3 || isempty(ord)
    ord = 4;
end

if nargin < 6 || isempty(verbose)
    verbose = true;
end

if ncomp > size(X0,2)
    ncomp = size(X0,2);
    warning('Number of components cannot exceed the dimensionality of X. Reducing ncomp to %i',ncomp);
end
dorand= nargin < 4 || isempty(a0);

maxiter = 500;
%%

m0 = mean(X0);
R0 = cov(X0);
[u,l,v] = svd(R0);
WM = u*sqrt(l)*v';
% WM = R0^(.5);
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

deflate = 1;
for dim = 1:ncomp

    R = X'*X;
    Rinv = pinv(R);


    if dorand
        a0 = deflate*randn(size(X0,2),1);
    end
   
    tol = 1e-6;
    d = Inf;

    a = WM*a0;

    % cm = cumulant(X*a,ord);

    iter = 1;
    clear ds kt sk
    if verbose
        nfp = fprintf('\nComponent %i, ',dim);
    end
    maxcum = 0;
    while d> tol && iter < maxiter

        r = X*a;


        if false && isnumeric(ord)
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
    if isnumeric(ord)
        cm = cumulant(X*anew,ord);
        if abs(cm)>abs(maxcum)
            iter_at_max = iter;
            maxcum=cm;
            maxa = anew;
        end
        if iter-iter_at_max  > 100
            iter = maxiter;
            fprintf('Max cumulant hasn''t been attained in 100 iterations. Stopping early.')
        end
    end
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
    a = sign(nansum(nl(X*anew)))*a; %Enforce positivity of the objective function extremum.
    A(:,dim) = a;   
    %%
    deflate = deflate*(eye(size(X,2))-a*a');
     X = X*deflate;
     if verbose
         switch ord
            case 3
                 nfp = fprintf('  iter %i, final skewness: %0.1f',iter,cumulant(r,ord));
            case 4
                nfp = fprintf('  iter %i, final ex. kurtosis: %0.1f',iter,cumulant(r,ord));
            otherwise
                nfp = fprintf('  iter %i, final ex. %ith order cumulant: %0.1f',iter,ord,cumulant(r,ord));
       
         end
     end
            
end

UM = pinv(WM)*A;
UM = UM.*sign(mean(dnl(X0*UM)));
