function out = two_gamma_mixture_em(z, varargin)
%TWO_GAMMA_MIXTURE_EM  EM for a two-component gamma mixture with shared shape.
%
% Model (motivated by the Gaussian state-space derivation of HOSD's
% smoothed-MS statistic X(t)^2):
%
%   z = X(t)^2 ~ pi_s * Gamma(alpha, beta_s) + (1 - pi_s) * Gamma(alpha, beta_n)
%
% Under x(t)/sigma_loc(t) ~ N(0, c_z^2) where z in {0,1} indexes state,
% and X(t)^2 a kernel-weighted sum of standardized squared samples, the
% shape parameter alpha is the effective DF of the smoothing kernel and
% the scales beta_n, beta_s scale with the state-conditional variance.
%
% INPUTS
%   z        N x 1 non-negative samples (= xrsm.^2 in the HOSD pipeline)
%   Optional name-value pairs:
%     'sample_w'    N x 1 non-negative weights (default ones; used by the
%                   caller to down-weight kernel-contaminated samples
%                   exactly as in detection_stats's frozen-noise fit).
%     'max_iter'    EM iterations (default 200)
%     'tol_logL'    convergence tolerance on log-likelihood (default 1e-7)
%     'init'        struct with .alpha, .beta_n, .beta_s, .pi_s
%     'outlier_log_pdf'  N x 1 vector of log f_o(z_t) for a FIXED outlier
%                   component (e.g. a frozen GPD log-density). When
%                   provided, EM runs as a 3-component mixture
%                       pi_n * f_n + pi_s * f_s + pi_o * f_o
%                   with f_o frozen and pi_o estimated subject to a cap
%                   (see outlier_pi_max). Outliers do not contribute to
%                   the noise / signal parameter updates -- the noise
%                   and signal beta_*/alpha estimators use only the
%                   non-outlier responsibilities. Default [] -> 2-component.
%     'outlier_pi_max'   cap on the outlier prior pi_o (default 0.01).
%                   The unconstrained EM estimate of pi_o is clipped at
%                   this value; the remaining mass (1 - pi_o) is
%                   distributed between pi_n and pi_s in proportion to
%                   their unconstrained responsibility-weighted counts.
%
% OUTPUTS
%   out.alpha        shared shape
%   out.beta_n       noise scale (smaller component)
%   out.beta_s       signal scale (larger component)
%   out.pi_s         signal prior probability
%   out.pi_o         outlier prior probability (0 when no outlier component)
%   out.gamma_t      N x 1 per-sample posterior P(signal | z_t)
%   out.gamma_o_t    N x 1 per-sample posterior P(outlier | z_t)
%   out.logL_path    per-iteration log-likelihood
%   out.n_iter       number of iterations performed
%   out.converged    logical

p = inputParser;
p.addParameter('sample_w', [], @isnumeric);
p.addParameter('max_iter', 200, @isnumeric);
p.addParameter('tol_logL', 1e-7, @isnumeric);
p.addParameter('init', [], @(s) isempty(s) || isstruct(s));
p.addParameter('outlier_log_pdf', [], @isnumeric);
p.addParameter('outlier_pi_max',  0.01, @isnumeric);
p.parse(varargin{:});
opt = p.Results;

z = z(:);
keep = z > 0 & isfinite(z);
z = z(keep);
N = numel(z);
if isempty(opt.sample_w)
    w = ones(N, 1);
else
    w = opt.sample_w(:);
    w = w(keep);
end
W = sum(w);

use_outlier = ~isempty(opt.outlier_log_pdf);
if use_outlier
    log_f_o = opt.outlier_log_pdf(:);
    assert(numel(log_f_o) >= numel(keep), ...
        'outlier_log_pdf must be at least as long as the unfiltered z input.');
    log_f_o = log_f_o(keep);
    pi_o = min(0.005, opt.outlier_pi_max);
else
    log_f_o = [];
    pi_o = 0;
end

% Initialization
if isempty(opt.init)
    % weighted method-of-moments on z to get a starting point;
    % then split the mode by a quantile-based heuristic
    m1   = sum(w .* z) / W;
    m2   = sum(w .* z.^2) / W;
    varz = max(m2 - m1^2, eps);
    a0   = max(m1^2 / varz, 0.5);
    b0   = varz / m1;
    pi_s = 0.05;
    % heuristic: split scale by 4x for signal so EM has clear basin of attraction
    alpha  = a0;
    beta_n = b0;
    beta_s = b0 * 4;
else
    alpha  = opt.init.alpha;
    beta_n = opt.init.beta_n;
    beta_s = opt.init.beta_s;
    pi_s   = opt.init.pi_s;
end

logz = log(z);
logL_path = zeros(opt.max_iter, 1);
converged = false;

for it = 1:opt.max_iter
    % --- E step: per-sample log-densities, log-sum-exp for stability ---
    log_f_n  = (alpha - 1) * logz - z / beta_n - alpha * log(beta_n) - gammaln(alpha);
    log_f_s  = (alpha - 1) * logz - z / beta_s - alpha * log(beta_s) - gammaln(alpha);
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

    % --- M step ---
    W_s = sum(w .* gamma_t);
    W_o = sum(w .* gamma_o_t);
    W_n = W - W_s - W_o;
    if use_outlier
        pi_o_unc = W_o / W;
        if pi_o_unc <= opt.outlier_pi_max
            pi_o_new = pi_o_unc;
            pi_s_new = W_s / W;
        else
            % Cap binds: pi_o pinned at outlier_pi_max; redistribute the
            % rest between n and s in their unconstrained-W ratio.
            pi_o_new = opt.outlier_pi_max;
            W_ns = max(W_n + W_s, eps);
            pi_s_new = (W_s / W_ns) * (1 - pi_o_new);
        end
    else
        pi_o_new = 0;
        pi_s_new = W_s / W;
    end

    % Component-conditional weighted means of z (responsibility-weighted).
    % Note: the alpha / beta updates use only the noise+signal mass; the
    % outlier component is held fixed and does not enter the M-step for
    % the gamma parameters.
    M_s = sum(w .* gamma_t          .* z) / max(W_s, eps);
    M_n = sum(w .* (1 - gamma_t - gamma_o_t) .* z) / max(W_n, eps);
    % Joint M-step optimality for alpha with beta_k = M_k/alpha:
    %   psi(alpha) - log(alpha) =
    %      [sum_t w_t (gamma_n + gamma_s) log z_t] / (W_n + W_s)
    %        - (W_n / (W_n+W_s)) log M_n - (W_s / (W_n+W_s)) log M_s
    W_ns = max(W_n + W_s, eps);
    rhs = sum(w .* (1 - gamma_o_t) .* logz) / W_ns ...
        - (W_n / W_ns) * log(max(M_n, eps)) ...
        - (W_s / W_ns) * log(max(M_s, eps));
    f_alpha = @(a) psi(a) - log(a) - rhs;
    try
        alpha_new = fzero(f_alpha, [1e-3, 1e4]);
        if isnan(alpha_new) || ~isfinite(alpha_new); alpha_new = alpha; end
    catch
        alpha_new = alpha;
    end
    beta_n_new = M_n / alpha_new;
    beta_s_new = M_s / alpha_new;

    % Ensure ordering beta_n < beta_s (noise has smaller scale).  If
    % they swap mid-iteration, relabel and re-run E next iter naturally.
    if beta_s_new < beta_n_new
        [beta_n_new, beta_s_new] = deal(beta_s_new, beta_n_new);
        pi_s_new = max(1 - pi_s_new - pi_o_new, 0);
    end

    % Convergence check on log-likelihood
    if it > 1 && abs(logL_path(it) - logL_path(it-1)) < opt.tol_logL * max(abs(logL_path(it)), eps)
        converged = true;
        alpha  = alpha_new;
        beta_n = beta_n_new;
        beta_s = beta_s_new;
        pi_s   = pi_s_new;
        pi_o   = pi_o_new;
        break
    end
    alpha  = alpha_new;
    beta_n = beta_n_new;
    beta_s = beta_s_new;
    pi_s   = pi_s_new;
    pi_o   = pi_o_new;
end

% Final E step for the reported posterior
log_f_n = (alpha - 1) * logz - z / beta_n - alpha * log(beta_n) - gammaln(alpha);
log_f_s = (alpha - 1) * logz - z / beta_s - alpha * log(beta_s) - gammaln(alpha);
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

% Pack output
out = struct();
out.alpha       = alpha;
out.beta_n      = beta_n;
out.beta_s      = beta_s;
out.pi_s        = pi_s;
out.pi_o        = pi_o;
out.gamma_t     = gamma_t;
out.gamma_o_t   = gamma_o_t;
out.logL_path   = logL_path(1:it);
out.n_iter      = it;
out.converged   = converged;
% SNR-equivalent ratio (squared, on z = X^2 scale)
out.snr2_ratio  = beta_s / beta_n;
% Noise/signal modes in x = sqrt(z) space (where x ~ X(t))
out.mode_x_n    = sqrt(max(beta_n * (2*alpha - 1) / 2, 0));
out.mode_x_s    = sqrt(max(beta_s * (2*alpha - 1) / 2, 0));
end
