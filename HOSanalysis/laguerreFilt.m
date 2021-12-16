function [Xfilt,structout] =laguerreFilt(X,ord,tau,fs)

% [Xfilt,str] =laguerreFilt(X,ord,alpha)
%Filter with laguerre polynomials up to specified order.
%
% Inputs:
%
%    X     : Input signal
%    ord   : Polynomial order as scalar value.
%    tau   : Time decay parameter. 
% 
%    Alternatively, the 2nd argument can be a struct array with fields
%               .order ->  polynomial order
%               .tau   ->  time decay parameter
% 
%    If X is a scalar integer, then the impulse response to the given 
%    number of samples is returned.
%    
%    
%    
% Outputs:
%
%    Xfilt : Filtered X 
%    str   : Struct array with order, tau and filter coeffificients.
%    
% C. Kovach 2019

if nargin < 4 || isempty(fs)    
    fs = 1;
end

if nargin < 3 || isempty(tau)
    tau = 1./log(2)./fs;
end

if nargin < 2
    ord = 0;
end

if isstruct(ord)
    structin = ord;
    if nargin>2 && ~isfield(structin,'fs')
       for k = 1:length(structin)
          structin(k).fs = tau; 
       end
    end
        
else
    structin = struct('order',ord,'tau',tau,'fs',fs);
end

if isscalar(X) % If the input is a scalar value then return the impulse response to specified order.
    X = [1;zeros(X-1,1)];
end

Xfilt = zeros(size(X,1),sum([structin.order]+1));

coli = 0;
for k = 1:length(structin)
    str = structin(k);
    
    alpha = exp(-1./str.tau/str.fs);

    B = [1 -1/alpha];
    A = [1 -alpha];
%     clear filtcoef
    coli = coli+1;
    Xfilt(:,coli) = filter([0 sqrt(1-alpha.^2)],[1 -alpha],full(X));
    for i = 1:str.order
        

        coli = coli+1;

        Xfilt(:,coli) = filter(B,A,Xfilt(:,coli-1));
        
      
%         filtcoef(i+1).order=i;
        
      
    end
    
%     str.filtcoefs=filtcoef;
    if nargout > 1
      str.A = A;
      str.B = B;
      structout(k) = str; %#ok<*AGROW>
    end
end


