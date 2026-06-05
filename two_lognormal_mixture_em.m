function out = two_lognormal_mixture_em(x, varargin)
%TWO_LOGNORMAL_MIXTURE_EM  Joint EM for noise + signal both lognormal on x.
%
% Model on the smoothed-MS statistic xrsm = X(t):
%   noise:  x ~ LogN(mu_n, sigma_n)
%   signal: x ~ LogN(mu_s, sigma_s)
%   x ~ pi_s * f_s(x) + (1 - pi_s - pi_o) * f_n(x) [ + pi_o * f_o(x) ]
%
% The frozen-noise variant (signal-only lognormal EM with noise pinned by
% the kernel-weighted MLE on sample_w) lives inline in detection_stats.m
% under the freeze_noise=true + signal_dist='lognormal' branch. This file
% provides the unfrozen / joint variant where both lognormals are updated
% each iteration -- analogous to TWO_GAMMA_MIXTURE_EM but in (mu, sigma)
% parameter space.
%
% Identifiability: lognormals are exchangeable, so we enforce mu_n < mu_s
% at the end of every M-step. If the M-step crosses (signal mu drops
% below noise mu), we swap component labels and rewrite pi_s -> 1 - pi_s
% - pi_o (same trick TWO_GAMMA_MIXTURE_EM uses on beta). In practice the
% asymmetric initialisation (mu_n at the 30th log-quantile, mu_s at the
% 70th) gives the right basin of attraction on every MASS subject tested
% with xrsm; the swap is a safety net for adversarial inputs.
%
% INPUTS
%   x        N x 1 non-negative samples (xrsm in HOSD pipeline). Samples
%            below an internal floor (1e-6) are dropped. In the standard
%            DETECTION_STATS dispatch this is pre-binned (log-spaced
%            x-bins, weights = bin counts via 'sample_w') so the EM
%            iterates over ~K bin centers instead of millions of raw
%            samples -- mathematically equivalent up to discretization
%            (see bin_em_K in detection_stats; ~0.03%-level parameter
%            agreement at K=500 log bins on heavy-tailed xrsm). The
%            same convention applies to every sibling helper
%            (TWO_GAMMA_MIXTURE_EM, GAMMA_LOGNORMAL_MIXTURE_EM,
%            GAMMA_NONCENTRAL_MIXTURE_EM).
%   Optional name-value pairs:
%     'sample_w'    N x 1 non-negative sample weights (default ones).
%                   In the binned-EM convention above this carries the
%                   bin counts -- so for two log-spaced bin centers x_i
%                   the weighted-MLE M-step recovers the same mu, sigma
%                   the raw EM would find on the unbinned samples.
%     'max_iter'    EM iterations (default 200)
%     'tol_logL'    relative tol on log-likelihood (default 1e-7)
%     'init'        struct with .mu_n, .sigma_n, .mu_s, .sigma_s, .pi_s
%     'outlier_log_pdf'  N x 1 frozen outlier log-pdf in x-space (e.g.
%                   GPD); when provided, EM runs as a 3-component mixture
%                   with the cap outlier_pi_max applied to pi_o (same
%                   semantics as the GAMMA_LOGNORMAL_MIXTURE_EM and
%                   TWO_GAMMA_MIXTURE_EM siblings).
%     'outlier_pi_max'   cap on pi_o when outlier_log_pdf is provided
%                        (default 0.01).
%
% OUTPUTS
%   out.mu_n, .sigma_n        noise lognormal params on x
%   out.mu_s, .sigma_s        signal lognormal params on x
%   out.pi_s, .pi_o           component priors (pi_o = 0 when 2-component)
%   out.gamma_t, .gamma_o_t   N x 1 per-sample posteriors P(signal|x), P(outlier|x)
%   out.logL_path             per-iteration log-likelihood
%   out.n_iter, .converged
%   out.mode_x_n, .mode_x_s   x-space modes (= exp(mu - sigma^2))

%C. Kovach 2026

p = inputParser;
p.addParameter('sample_w', [], @isnumeric);
p.addParameter('max_iter', 200, @isnumeric);
p.addParameter('tol_logL', 1e-7, @isnumeric);
p.addParameter('init', [], @(s) isempty(s) || isstruct(s));
p.addParameter('outlier_log_pdf', [], @isnumeric);
p.addParameter('outlier_pi_max',  0.01, @isnumeric);
p.parse(varargin{:});
opt = p.Results;

x = x(:);
keep = x > 1e-6 & isfinite(x);
x = x(keep);
N = numel(x);
if isempty(opt.sample_w)
    w = ones(N, 1);
else
    w = opt.sample_w(:); w = w(keep);
end
W = sum(w);

use_outlier = ~isempty(opt.outlier_log_pdf);
if use_outlier
    log_f_o = opt.outlier_log_pdf(:);
    assert(numel(log_f_o) >= numel(keep), ...
        'outlier_log_pdf must be at least as long as the unfiltered x input.');
    log_f_o = log_f_o(keep);
    pi_o = min(0.005, opt.outlier_pi_max);
else
    log_f_o = [];
    pi_o = 0;
end

logx = log(x);

% --- Initialization ---
if isempty(opt.init)
    % Asymmetric quantile-based seed: noise at lower log-x, signal at
    % upper. The 30/70 split (rather than 25/75) leaves a comfortable
    % gap that survives the first M-step when pi_s defaults to 0.10
    % (i.e. signal mass is small and the M-step nudge from the 70th
    % quantile is modest).
    q_lo = weighted_quantile(logx, w, 0.30);
    q_hi = weighted_quantile(logx, w, 0.70);
    mu_n = q_lo;
    mu_s = q_hi;
    % Half-window standard deviations of log x, floored at 0.3 to keep
    % both components from initialising too narrow.
    med = weighted_quantile(logx, w, 0.50);
    sel_lo = logx <= med;
    sel_hi = ~sel_lo;
    sigma_n = max(weighted_std(logx(sel_lo), w(sel_lo)), 0.3);
    sigma_s = max(weighted_std(logx(sel_hi), w(sel_hi)), 0.3);
    pi_s    = 0.10;
else
    mu_n    = opt.init.mu_n;
    sigma_n = opt.init.sigma_n;
    mu_s    = opt.init.mu_s;
    sigma_s = opt.init.sigma_s;
    pi_s    = opt.init.pi_s;
end

logL_path = zeros(opt.max_iter, 1);
converged = false;

for it = 1:opt.max_iter
    % --- E-step (x-space; Jacobian dz/dx is shared so cancels in responsibilities) ---
    log_f_n = -logx - log(sigma_n) - 0.5*log(2*pi) - (logx - mu_n).^2 / (2*sigma_n^2);
    log_f_s = -logx - log(sigma_s) - 0.5*log(2*pi) - (logx - mu_s).^2 / (2*sigma_s^2);

    log_pi_n = log(max(1 - pi_s - pi_o, eps));
    log_pi_s = log(max(pi_s, eps));
    log_w_n  = log_pi_n + log_f_n;
    log_w_s  = log_pi_s + log_f_s;
    if use_outlier
        log_pi_o = log(max(pi_o, eps));
        log_w_o  = log_pi_o + log_f_o;
        m_log = max(max(log_w_n, log_w_s), log_w_o);
        log_denom = m_log + log(exp(log_w_n - m_log) + exp(log_w_s - m_log) + exp(log_w_o - m_log));
        gamma_o_t = exp(log_w_o - log_denom);
    else
        m_log = max(log_w_n, log_w_s);
        log_denom = m_log + log(exp(log_w_n - m_log) + exp(log_w_s - m_log));
        gamma_o_t = zeros(N, 1);
    end
    gamma_t = exp(log_w_s - log_denom);
    logL_path(it) = sum(w .* log_denom);

    % --- M-step ---
    W_s = sum(w .* gamma_t);
    W_o = sum(w .* gamma_o_t);
    W_n = W - W_s - W_o;
    if use_outlier
        pi_o_unc = W_o / W;
        if pi_o_unc <= opt.outlier_pi_max
            pi_o_new = pi_o_unc;
            pi_s_new = W_s / W;
        else
            pi_o_new = opt.outlier_pi_max;
            W_ns = max(W_n + W_s, eps);
            pi_s_new = (W_s / W_ns) * (1 - pi_o_new);
        end
    else
        pi_o_new = 0;
        pi_s_new = W_s / W;
    end

    % Signal lognormal -- closed-form weighted MLE on log x with weight gamma_t.
    mu_s_new    = sum(w .* gamma_t .* logx) / max(W_s, eps);
    sigma_s_new = sqrt(max(sum(w .* gamma_t .* (logx - mu_s_new).^2) / max(W_s, eps), 1e-6));

    % Noise lognormal -- closed-form weighted MLE on log x with weight
    % (1 - gamma_t - gamma_o_t). Outliers do not contribute to noise updates.
    w_n     = w .* (1 - gamma_t - gamma_o_t);
    W_n_eff = max(sum(w_n), eps);
    mu_n_new    = sum(w_n .* logx) / W_n_eff;
    sigma_n_new = sqrt(max(sum(w_n .* (logx - mu_n_new).^2) / W_n_eff, 1e-6));

    % Enforce ordering mu_n < mu_s (noise at lower x-location). Mirror
    % of the beta_n < beta_s swap in TWO_GAMMA_MIXTURE_EM. If the M-step
    % crosses, swap labels (params + responsibilities) and rewrite
    % pi_s -> 1 - pi_s_new - pi_o_new. The next E-step then carries the
    % corrected labels forward.
    if mu_s_new < mu_n_new
        [mu_n_new, mu_s_new]       = deal(mu_s_new, mu_n_new);
        [sigma_n_new, sigma_s_new] = deal(sigma_s_new, sigma_n_new);
        pi_s_new = max(1 - pi_s_new - pi_o_new, 0);
    end

    % --- Convergence ---
    if it > 1 && abs(logL_path(it) - logL_path(it-1)) < opt.tol_logL * max(abs(logL_path(it)), eps)
        converged = true;
        mu_n    = mu_n_new;    sigma_n = sigma_n_new;
        mu_s    = mu_s_new;    sigma_s = sigma_s_new;
        pi_s    = pi_s_new;    pi_o    = pi_o_new;
        break
    end
    mu_n    = mu_n_new;    sigma_n = sigma_n_new;
    mu_s    = mu_s_new;    sigma_s = sigma_s_new;
    pi_s    = pi_s_new;    pi_o    = pi_o_new;
end

% --- Final E-step for the reported responsibilities ---
log_f_n = -logx - log(sigma_n) - 0.5*log(2*pi) - (logx - mu_n).^2 / (2*sigma_n^2);
log_f_s = -logx - log(sigma_s) - 0.5*log(2*pi) - (logx - mu_s).^2 / (2*sigma_s^2);
log_pi_n = log(max(1 - pi_s - pi_o, eps));
log_pi_s = log(max(pi_s, eps));
log_w_n  = log_pi_n + log_f_n;
log_w_s  = log_pi_s + log_f_s;
if use_outlier
    log_pi_o = log(max(pi_o, eps));
    log_w_o  = log_pi_o + log_f_o;
    m_log = max(max(log_w_n, log_w_s), log_w_o);
    log_denom = m_log + log(exp(log_w_n - m_log) + exp(log_w_s - m_log) + exp(log_w_o - m_log));
    gamma_o_t = exp(log_w_o - log_denom);
else
    m_log = max(log_w_n, log_w_s);
    log_denom = m_log + log(exp(log_w_n - m_log) + exp(log_w_s - m_log));
    gamma_o_t = zeros(N, 1);
end
gamma_t = exp(log_w_s - log_denom);

out = struct();
out.mu_n      = mu_n;
out.sigma_n   = sigma_n;
out.mu_s      = mu_s;
out.sigma_s   = sigma_s;
out.pi_s      = pi_s;
out.pi_o      = pi_o;
out.gamma_t   = gamma_t;
out.gamma_o_t = gamma_o_t;
out.logL_path = logL_path(1:it);
out.n_iter    = it;
out.converged = converged;
% Modes in x-space. Lognormal(mu, sigma) on x has mode exp(mu - sigma^2).
out.mode_x_n  = exp(mu_n - sigma_n^2);
out.mode_x_s  = exp(mu_s - sigma_s^2);
end

% =========================================================================
function q = weighted_quantile(x, w, p)
%WEIGHTED_QUANTILE  Linear-interp-free weighted quantile (step CDF).
[xs, ord] = sort(x(:));
ws = w(ord);
cw = cumsum(ws) / max(sum(ws), eps);
idx = find(cw >= p, 1, 'first');
if isempty(idx); idx = numel(xs); end
q = xs(idx);
end

% =========================================================================
function s = weighted_std(x, w)
%WEIGHTED_STD  Plain weighted std (no bias correction; used only for init).
W = max(sum(w), eps);
mu = sum(w .* x) / W;
s = sqrt(max(sum(w .* (x - mu).^2) / W, 0));
end
