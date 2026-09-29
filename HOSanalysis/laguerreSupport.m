function L = laguerreSupport(ord, tau, fs, nmax, tol)
%LAGUERRESUPPORT  Number of lags the discrete Laguerre functions occupy.
%
%   L = laguerreSupport(ord, tau, fs, nmax, tol)
%
%   The ord+1 Laguerre functions of laguerreFilt / laguerreBasis (time
%   constant tau in seconds, sampling rate fs in Hz) decay exponentially
%   but never vanish. L is the number of lags after which all of them stay
%   below tol times the largest magnitude of the set: lags 0 .. L-1 hold
%   everything that matters. The higher the order the longer the support
%   (roughly 4 (ord+1) tau before the exponential takes over).
%
%   nmax   cap on L (default Inf), e.g. the length of the record
%   tol    relative magnitude regarded as zero (default 1e-12)
%
%   Used to zero-pad the record so that the FFT convolution of an impulse
%   train with a Laguerre basis is LINEAR: an event near the end of the
%   record must not wrap its kernel onto the beginning (model.designMtx).
%
%   See also LAGUERREFILT, LAGUERREBASIS, MODEL.

% C. Kovach 2026

if nargin < 4 || isempty(nmax); nmax = Inf; end
if nargin < 5 || isempty(tol);  tol = 1e-12; end
if ~(isfinite(tau) && tau > 0 && isfinite(fs) && fs > 0 && isfinite(ord) && ord >= 0)
    error('laguerreSupport:badInput', 'ord must be >= 0, tau and fs positive and finite.');
end
Ltry = max(16, ceil(8 * (ord + 1) * tau * fs));
while true
    Ltry = min(Ltry, nmax);
    h = laguerreFilt(Ltry + 1, ord, tau, fs);   % row 1 is the filter's one-sample delay
    h = h(2:end, :);
    a = max(abs(h), [], 2);
    last = find(a > tol * max(a), 1, 'last');
    if isempty(last)
        L = 0;
        return
    end
    if last < Ltry || Ltry >= nmax
        L = min(last, nmax);
        return
    end
    Ltry = 2 * Ltry;
end
end
