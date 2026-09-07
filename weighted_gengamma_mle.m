function [a, d, p, info] = weighted_gengamma_mle(x, w, varargin)
%WEIGHTED_GENGAMMA_MLE  Generalized gamma (Stacy) MLE with sample weights.
%
%   [a, d, p, info] = weighted_gengamma_mle(x, w)
%
%   Fits the three-parameter generalized gamma density
%
%       f(x; a, d, p) = p / (a^d Gamma(d/p)) * x^(d-1) * exp(-(x/a)^p),  x > 0
%
%   (scale a > 0, shape d > 0, power p > 0) by weighted maximum
%   likelihood. The family nests the gamma on x (p = 1, shape d, scale
%   a), the Weibull (d = p), the gamma on x^2 used by detection_stats'
%   'gamma' family (p = 2, d = 2*shape) and, as limits, the lognormal
%   (d -> Inf) and the half-normal (p = 2, d = 1).
%
%   Weights are the usual binned-data convention (bin centres in x, bin
%   counts or noise responsibilities in w). For fixed (d, p) the scale is
%   closed-form,  a^p = p * sum(w x^p) / (d * sum(w)),  so the profile
%   likelihood is a 2-D problem in (log d, log p), solved with fminsearch
%   from two starts: the Weibull solution (d = p = k) and the gamma-on-x
%   solution (p = 1). Data are rescaled by their weighted geometric mean
%   so x^p cannot overflow.
%
%   info.negloglik   weighted negative log-likelihood at the optimum
%   info.start       'weibull' | 'gamma' (the start that won)
%   info.exitflag    fminsearch exit flag
%
%   Copyright (C) 2026 Christopher K. Kovach, University of Nebraska
%   Medical Center.

    keep = isfinite(x) & x > 0 & isfinite(w) & w > 0;
    x = double(x(keep)); w = double(w(keep));
    if numel(x) < 10
        error('weighted_gengamma_mle:tooFew', 'too few samples');
    end
    x = x(:); w = w(:);
    sw = sum(w);
    lx = log(x);
    ml = sum(w .* lx) / sw;
    c  = exp(ml);                       % weighted geometric mean
    y  = x / c; ly = lx - ml;           % mean_w(ly) = 0
    sd = sqrt(max(sum(w .* ly.^2) / sw, 1e-12));

    % Start 1: Weibull  (d = p = k, k from SD(log X) = pi / (k sqrt 6))
    k0 = 1.2825 / sd;
    th1 = [log(k0), log(k0)];
    % Start 2: gamma on x (p = 1), Minka's closed-form shape
    s = log(sum(w .* y) / sw) - 0;      % log(mean_w y) - mean_w(log y), the latter is 0
    if s > 1e-10
        d0 = (3 - s + sqrt((s - 3)^2 + 24 * s)) / (12 * s);
    else
        d0 = 1e3;
    end
    th2 = [log(d0), 0];

    opts = optimset('Display', 'off', 'TolX', 1e-7, 'TolFun', 1e-9, ...
                    'MaxFunEvals', 4000, 'MaxIter', 2000);
    if ~isempty(varargin); opts = optimset(opts, varargin{:}); end
    nll = @(th) negloglik(th, y, ly, w, sw);
    best = Inf; th = th1; startName = 'weibull'; flag = 0;
    for st = {{th1, 'weibull'}, {th2, 'gamma'}}
        try
            [thk, fk, ek] = fminsearch(nll, st{1}{1}, opts);
        catch
            continue
        end
        if isfinite(fk) && fk < best
            best = fk; th = thk; startName = st{1}{2}; flag = ek;
        end
    end
    d = exp(th(1)); p = exp(th(2));
    ay = (p * sum(w .* y.^p) / (d * sw))^(1 / p);
    a = c * ay;
    info = struct('negloglik', best, 'start', startName, 'exitflag', flag);
end

function v = negloglik(th, y, ly, w, sw)
    d = exp(th(1)); p = exp(th(2));
    if d < 0.02 || d > 1e5 || p < 0.02 || p > 100
        v = Inf; return
    end
    yp = y.^p;
    m  = sum(w .* yp) / sw;             % mean_w(y^p)
    if ~isfinite(m) || m <= 0; v = Inf; return; end
    ap = p * m / d;                     % profile a^p (in y units)
    % sum_i w_i log f(y_i) with a^p profiled out:
    %   sw * [log p - (d/p) log(a^p) - gammaln(d/p) - d/p] + (d - 1) sum(w ly)
    ll = sw * (log(p) - (d / p) * log(ap) - gammaln(d / p) - d / p) + (d - 1) * sum(w .* ly);
    v = -ll;
    if ~isfinite(v); v = Inf; end
end
