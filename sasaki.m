function out = sasaki(N)

%
% Minimum bias window for bispectral estimation following Sasaki, Sato and
% Yamashita (1975).
%
% Sasaki, K., T. Sato, and Y. Yamashita (1975). Minimum bias windows for 
%   bispectral estimation. Journal of Sound and Vibration 40(1), 139?148.
%


t = linspace(-1,1,N)';

saswin = @(x)1/pi.*abs(sin(pi*x)) + (1-abs(x)).*cos(pi*x);

out = saswin(t);

tol = 1e-9;

out = round(out./tol)*tol;