function out = gamma_lognormal_mixture_em(x, varargin)
%GAMMA_LOGNORMAL_MIXTURE_EM  EM for noise (gamma-on-z) + signal (lognormal-on-x).
%
% Model on the smoothed-MS statistic xrsm = X(t):
%   noise:  z = x^2 ~ Gamma(alpha, beta_n)  -> f_n(x) = 2x * gampdf(x^2; alpha, beta_n)
%   signal: x       ~ LogN(mu_s, sigma_s)   -> f_s(x) = lognpdf(x; mu_s, sigma_s)
%   x ~ pi_s * f_s(x) + (1-pi_s) * f_n(x)
%
% The lognormal signal family is sub-exponential -- its right tail decays as
% exp(-(log x - mu)^2 / 2 sigma^2), much slower than gamma's exp(-x^2/beta_n).
% This was added after the 2-gamma joint EM was found to be too light-tailed
% to fit empirical xrsm distributions (a GPD outlier component was eating
% real-spindle mass; see the relabeling-test analysis).
%
% INPUTS
%   x        N x 1 non-negative samples (xrsm in HOSD pipeline). Samples
%            below an internal floor (1e-6) are dropped.
%   Optional name-value pairs:
%     'sample_w'    N x 1 non-negative sample weights (default ones)
%     'max_iter'    EM iterations (default 200)
%     'tol_logL'    relative tol on log-likelihood (default 1e-7)
%     'init'        struct with .alpha, .beta_n, .mu_s, .sigma_s, .pi_s
%     'outlier_log_pdf'  N x 1 frozen outlier log-pdf (e.g. GPD); when
%                   provided, EM runs as a 3-component mixture with the
%                   cap outlier_pi_max applied to pi_o (same semantics as
%                   two_gamma_mixture_em). Default [] -> 2-component.
%     'outlier_pi_max'   cap on pi_o when outlier_log_pdf is provided
%                        (default 0.01).
%
% OUTPUTS
%   out.alpha, .beta_n        gamma noise params on z = x^2
%   out.mu_s, .sigma_s        lognormal signal params on x
%   out.pi_s, .pi_o           component priors (pi_o = 0 when 2-component)
%   out.gamma_t, .gamma_o_t   N x 1 per-sample posteriors P(signal|x), P(outlier|x)
%   out.logL_path             per-iteration log-likelihood
%   out.n_iter, .converged

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
    log_f_o = log_f_o(keep);
    pi_o = min(0.005, opt.outlier_pi_max);
else
    log_f_o = [];
    pi_o = 0;
end

logx = log(x);
z    = x.^2;
logz = log(z);

% Initialisation
if isempty(opt.init)
    % Method-of-moments on z for noise; lognormal init from log-quantile spread
    m1 = sum(w .* z) / W; m2 = sum(w .* z.^2) / W;
    varz = max(m2 - m1^2, eps);
    alpha  = max(m1^2/varz, 0.5);
    beta_n = m1 / max(alpha, eps);   % noise covers the bulk
    % Signal init: log-mean and log-sd from the upper-half of samples
    mu_s    = quantile(logx, 0.6);
    sigma_s = max(std(logx(logx >= quantile(logx, 0.5))), 0.3);
    pi_s    = 0.10;
else
    alpha   = opt.init.alpha;
    beta_n  = opt.init.beta_n;
    mu_s    = opt.init.mu_s;
    sigma_s = opt.init.sigma_s;
    pi_s    = opt.init.pi_s;
end

logL_path = zeros(opt.max_iter, 1);
converged = false;

for it = 1:opt.max_iter
    % --- E-step ---
    % noise log-pdf in x-space: f_n(x) = 2x * gampdf(x^2; alpha, beta_n)
    log_f_n = log(2) + logx + (alpha-1)*logz - z/beta_n - alpha*log(beta_n) - gammaln(alpha);
    % signal log-pdf in x-space: lognormal
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

    % Signal lognormal -- closed-form weighted MLE on log x with weight gamma_t
    mu_s_new    = sum(w .* gamma_t .* logx) / max(W_s, eps);
    sigma_s_new = sqrt(max(sum(w .* gamma_t .* (logx - mu_s_new).^2) / max(W_s, eps), 1e-6));

    % Noise gamma on z with weight (1 - gamma_t - gamma_o_t).
    % Weighted MoM-init for alpha, then Newton on log(a) - psi(a) = s_n.
    w_n   = w .* (1 - gamma_t - gamma_o_t);
    sum_wn = max(sum(w_n), eps);
    mean_z_n   = sum(w_n .* z)   / sum_wn;
    mean_logz_n = sum(w_n .* logz) / sum_wn;
    s_n = log(max(mean_z_n, eps)) - mean_logz_n;       % >= 0; 0 means delta-distribution
    % closed-form approximation for alpha (Minka 2002), then refine
    a_init = (3 - s_n + sqrt((s_n - 3)^2 + 24*s_n)) / max(12*s_n, eps);
    a = max(a_init, 1e-3);
    for ni = 1:50
        f  = log(a) - psi(a) - s_n;
        fp = 1/a - psi(1, a);
        a_new = a - f/fp;
        if ~isfinite(a_new) || a_new <= 0, break, end
        if abs(a_new - a) < 1e-7 * abs(a) + 1e-10, a = a_new; break, end
        a = a_new;
    end
    alpha_new  = a;
    beta_n_new = mean_z_n / max(alpha_new, eps);

    % Convergence
    if it > 1 && abs(logL_path(it) - logL_path(it-1)) < opt.tol_logL * max(abs(logL_path(it)), eps)
        converged = true;
        alpha   = alpha_new;   beta_n  = beta_n_new;
        mu_s    = mu_s_new;    sigma_s = sigma_s_new;
        pi_s    = pi_s_new;    pi_o    = pi_o_new;
        break
    end
    alpha   = alpha_new;   beta_n  = beta_n_new;
    mu_s    = mu_s_new;    sigma_s = sigma_s_new;
    pi_s    = pi_s_new;    pi_o    = pi_o_new;
end

% Final E-step for the reported responsibilities
log_f_n = log(2) + logx + (alpha-1)*logz - z/beta_n - alpha*log(beta_n) - gammaln(alpha);
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
out.alpha     = alpha;
out.beta_n    = beta_n;
out.mu_s      = mu_s;
out.sigma_s   = sigma_s;
out.pi_s      = pi_s;
out.pi_o      = pi_o;
out.gamma_t   = gamma_t;
out.gamma_o_t = gamma_o_t;
out.logL_path = logL_path(1:it);
out.n_iter    = it;
out.converged = converged;
% Noise mode in x-space (z = x^2): x_mode = sqrt(beta_n (2 alpha - 1)/2) for alpha > 1/2
out.mode_x_n = sqrt(max(beta_n * (2*alpha - 1) / 2, 0));
% Lognormal mode in x-space: exp(mu - sigma^2)
out.mode_x_s = exp(mu_s - sigma_s^2);
end
