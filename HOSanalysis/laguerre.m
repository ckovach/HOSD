function [L,A] =laguerre(nx,ord,alpha)

%Laguerre polynomial

if nargin < 3
    alpha = 1;
end

if nargin < 2
    ord = 0;
end

L = zeros(nx,length(ord));
F = [1;zeros(nx-1,1)];

rootsA = sqrt(alpha);
rootsB = [];

for i = 0:ord
    A = poly(rootsA);
    B = sqrt(1-alpha)*poly(rootsB);
    L(:,i+1) = filter(B,A,F);
    rootsA(end+1)= sqrt(alpha);
    rootsB(end+1) = 1/sqrt(alpha);
end
   


