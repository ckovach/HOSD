function out = gamma_noncentral_mixture_em(z, varargin)
%GAMMA_NONCENTRAL_MIXTURE_EM  EM for noise (central gamma) + signal (non-central gamma).
%
% Model:
%   noise:  z | state=0 ~ Gamma(alpha, beta_n)
%   signal: z | state=1 ~ NonCentralGamma(alpha, beta_s, lambda)
%   z ~ pi_s * NCgamma(alpha, beta_s, lambda) + (1-pi_s) * Gamma(alpha, beta_n)
%
% Implementation uses the Poisson-augmentation form of NCgamma:
%   NCgamma(z; alpha, beta, lambda) = exp(-lambda/2) * sum_k Poisson(k; lambda/2) *
%                                       Gamma_pdf(z; alpha + k, beta).
% Augmenting with the latent count K | Z ~ r_k(z) yields closed-form M-step
% updates for pi_s, beta_n, beta_s, lambda; alpha is solved by a 1D numerical
% root-find using the proper non-central digamma equation
%   (W_n/W) psi(alpha) + (1/W) sum_t gamma_t E[psi(alpha+K_t) | z_t]
%       = mean_w(log z) - (W_n/W) log beta_n - (W_s/W) log beta_s
% with the augmentation posterior r_k(z_t) computed from the previous-iter
% parameters (standard EM-with-augmentation pattern).
%
% INPUTS
%   z         N x 1 non-negative samples
%   Optional name-value:
%     'mode'        'fixed_scale' -> beta_n = beta_s = beta (only lambda
%                                    differentiates signal from noise);
%                                    closed-form combined-beta M-step.
%                   'free_scale'  -> beta_n, beta_s, lambda all free (default)
%     'max_iter'    EM iterations (default 200)
%     'tol_logL'    relative convergence tol on log-likelihood (default 1e-6)
%     'init'        struct with .alpha, .beta_n, .beta_s, .lambda, .pi_s
%     'sample_w'    per-sample weights (default ones)
%     'K_max'       truncation of the Poisson-augmentation sum (default
%                   chosen per-iteration as ceil(lambda/2 + 4*sqrt(lambda/2+1)))
%     'verbose'     true to print per-iteration parameters
%
% OUTPUTS
%   out.mode, .alpha, .beta_n, .beta_s, .lambda, .pi_s
%   out.gamma_t       N x 1 P(signal | z) at convergence
%   out.logL_path     per-iter log-likelihood
%   out.n_iter, .converged

p = inputParser;
p.addParameter('mode', 'free_scale', @ischar);
p.addParameter('max_iter', 200, @isnumeric);
p.addParameter('tol_logL', 1e-6, @isnumeric);
p.addParameter('init', [], @(s) isempty(s) || isstruct(s));
p.addParameter('sample_w', [], @isnumeric);
p.addParameter('K_max', [], @isnumeric);
p.addParameter('verbose', false, @islogical);
p.parse(varargin{:});
opt = p.Results;

assert(any(strcmp(opt.mode, {'fixed_scale', 'free_scale'})), ...
    'mode must be ''fixed_scale'' or ''free_scale''');

z = z(:);
keep = z > 0 & isfinite(z);
z = z(keep);
N = numel(z);
if isempty(opt.sample_w)
    w = ones(N, 1);
else
    w = opt.sample_w(:); w = w(keep);
end
W = sum(w);
logz = log(z);

% Initialization
if isempty(opt.init)
    m1 = sum(w .* z) / W;
    m2 = sum(w .* z.^2) / W;
    varz = max(m2 - m1^2, eps);
    alpha = max(m1^2 / varz, 0.5);
    beta_n = varz / m1;
    beta_s = beta_n;
    lambda = 4 * alpha;
    pi_s = 0.10;
else
    alpha  = opt.init.alpha;
    beta_n = opt.init.beta_n;
    beta_s = opt.init.beta_s;
    lambda = opt.init.lambda;
    pi_s   = opt.init.pi_s;
end
if strcmp(opt.mode, 'fixed_scale'); beta_s = beta_n; end

logL_path = zeros(opt.max_iter, 1);
converged = false;

for it = 1:opt.max_iter
    % --- E-step: noise pdf (closed form) ---
    log_f_n = (alpha - 1) * logz - z/beta_n - alpha*log(beta_n) - gammaln(alpha);

    % --- E-step: signal pdf + Poisson-augmentation posterior ---
    [E_K, E_psi_aK, r_k, log_f_s, K_vals] = ...
        noncentral_aug_posterior(z, alpha, beta_s, lambda, opt.K_max);

    % Mixture responsibility for state (signal vs noise)
    log_pi_n = log(max(1 - pi_s, eps));
    log_pi_s = log(max(pi_s,     eps));
    m_log    = max(log_pi_n + log_f_n, log_pi_s + log_f_s);
    log_denom = m_log + log(exp(log_pi_n + log_f_n - m_log) ...
                          + exp(log_pi_s + log_f_s - m_log));
    gamma_t  = exp(log_pi_s + log_f_s - log_denom);
    logL_path(it) = sum(w .* log_denom);

    % --- M-step ---
    W_s = sum(w .* gamma_t);
    W_n = W - W_s;
    pi_s_new = W_s / W;

    % Sufficient statistics for the augmentation-weighted updates
    sum_gz   = sum(w .* gamma_t .* z);
    sum_gEK  = sum(w .* gamma_t .* E_K);
    sum_gaEK = sum(w .* gamma_t .* (alpha + E_K));

    % Lambda: closed-form. lambda_new = 2 * (responsibility-weighted mean
    % of E[K|z]) in the signal component.
    lambda_new = 2 * sum_gEK / max(W_s, eps);

    if strcmp(opt.mode, 'fixed_scale')
        % Shared beta_n = beta_s = beta. Combined M-step:
        %   d/db [Σ (1-γ) log_gam(z;α,b) + Σ γ log_NCgam(z;α,b,λ)] = 0
        %   => β = Σ w z / [α W + Σ w γ E[K|z]]
        beta_shared = sum(w .* z) / max(alpha * W + sum_gEK, eps);
        beta_n_new = beta_shared;
        beta_s_new = beta_shared;
    else
        % beta_n: classic central-gamma weighted MLE conditional on alpha
        beta_n_new = sum(w .* (1 - gamma_t) .* z) / max(W_n * alpha, eps);
        % beta_s: closed-form for non-central via augmentation
        %   β_s = Σ γ z / Σ γ (α + E[K|z])
        beta_s_new = sum_gz / max(sum_gaEK, eps);
    end

    % alpha: numerical root-find on the non-central digamma equation
    %   (W_n/W) ψ(α) + (W_s/W) E_avg[ψ(α + K)] = mean_w(log z)
    %                                            - (W_n/W) log β_n_new
    %                                            - (W_s/W) log β_s_new
    % r_k(z_t) is treated as fixed at the OLD parameters (standard EM-
    % augmentation pattern); only ψ(α + k) varies with the new α.
    v_k = (w .* gamma_t).' * r_k;       % 1 x (K_max+1), weighted r_k sum over signal
    rhs = (sum(w .* logz) / W) ...
        - (W_n / W) * log(max(beta_n_new, eps)) ...
        - (W_s / W) * log(max(beta_s_new, eps));
    K_vals_row = K_vals(:).';
    f_alpha = @(a) (W_n/W) * psi(a) + (1/W) * sum(v_k .* psi(a + K_vals_row)) - rhs;
    try
        alpha_new = fzero(f_alpha, [1e-3, 1e4]);
        if ~isfinite(alpha_new); alpha_new = alpha; end
    catch
        alpha_new = alpha;
    end

    if opt.verbose
        fprintf('  it %3d  logL=%.4e  alpha=%.3f  beta_n=%.3f  beta_s=%.3f  lambda=%.3f  pi_s=%.4f\n', ...
            it, logL_path(it), alpha_new, beta_n_new, beta_s_new, lambda_new, pi_s_new);
    end

    if it > 1 && abs(logL_path(it) - logL_path(it-1)) < opt.tol_logL * max(abs(logL_path(it)), eps)
        converged = true;
        alpha  = alpha_new;
        beta_n = beta_n_new; beta_s = beta_s_new;
        lambda = lambda_new; pi_s = pi_s_new;
        break
    end
    alpha  = alpha_new;
    beta_n = beta_n_new; beta_s = beta_s_new;
    lambda = lambda_new; pi_s = pi_s_new;
end

% Final E-step for the reported posterior
log_f_n = (alpha - 1) * logz - z/beta_n - alpha*log(beta_n) - gammaln(alpha);
[~, ~, ~, log_f_s, ~] = noncentral_aug_posterior(z, alpha, beta_s, lambda, opt.K_max);
log_pi_n = log(max(1-pi_s, eps)); log_pi_s = log(max(pi_s, eps));
m_log = max(log_pi_n + log_f_n, log_pi_s + log_f_s);
log_denom = m_log + log(exp(log_pi_n + log_f_n - m_log) ...
                      + exp(log_pi_s + log_f_s - m_log));
gamma_t = exp(log_pi_s + log_f_s - log_denom);

out = struct( ...
    'mode',      opt.mode, ...
    'alpha',     alpha, ...
    'beta_n',    beta_n, ...
    'beta_s',    beta_s, ...
    'lambda',    lambda, ...
    'pi_s',      pi_s, ...
    'gamma_t',   gamma_t, ...
    'logL_path', logL_path(1:it), ...
    'n_iter',    it, ...
    'converged', converged);
end

% ============================================================================
function [E_K, E_psi_aK, r_k, log_ncgam, K_vals] = ...
    noncentral_aug_posterior(z, alpha, beta, lambda, K_max)
%NONCENTRAL_AUG_POSTERIOR  Poisson-augmentation posterior for NCgamma.
%
% Returns per-sample:
%   E_K        N x 1   E[K | z]
%   E_psi_aK   N x 1   E[ψ(α + K) | z]
%   r_k        N x (K_max+1)  posterior P(K=k | z) — needed by the alpha M-step
%   log_ncgam  N x 1   log f_NCgam(z; α, β, λ)
%   K_vals     1 x (K_max+1)  the integer K grid used
%
% NCgamma pdf in Poisson-mixture form (note exp(-lambda/2) pulled out):
%   f_NCgam(z) = exp(-λ/2) Σ_{k>=0} (λ/2)^k / k! · f_gam(z; α + k, β)
% log_a_k(z) := log[(λ/2)^k / k! · f_gam(z; α+k, β)]
%             = k log(λ/2) - log k! + (α+k-1) log z - z/β - (α+k) log β - log Γ(α+k)
% Then log f_NCgam = -λ/2 + log Σ_k exp(log_a_k), and
%      r_k(z) = exp(log_a_k(z) - log Σ_k exp(log_a_k(z))).
% Numerically stable via log-sum-exp on log_a_k.

if nargin < 5 || isempty(K_max)
    K_max = max(5, ceil(lambda/2 + 4 * sqrt(lambda/2 + 1)));
end
z = z(:);
N = numel(z);
K_vals = 0:K_max;

log_z = log(max(z, eps));
log_b = log(max(beta, eps));
log_lh = log(max(lambda/2, eps));

% k-only terms: log(λ/2)^k / k! - log Γ(α+k)
k_only = K_vals * log_lh - gammaln(K_vals + 1) - gammaln(alpha + K_vals);   % 1 x (K+1)

% mixed terms: k * (log z - log β); broadcast over t and k
%   bracket(t, k+1) = k_only(k+1) + k * (log z(t) - log β)
bracket = K_vals .* (log_z - log_b);             % N x (K+1) (implicit expansion)
bracket = bracket + k_only;                      % N x (K+1) (implicit expansion)

% t-only common: (α-1) log z(t) - z(t)/β - α log β
common = (alpha - 1) * log_z - z/beta - alpha*log_b;

log_a = bracket + common;                        % N x (K+1)

% Stable normalisation per row
m = max(log_a, [], 2);
ln = log_a - m;
a  = exp(ln);
sum_a = sum(a, 2);
r_k = a ./ sum_a;

E_K       = r_k * K_vals(:);
psi_vals  = psi(alpha + K_vals);
E_psi_aK  = r_k * psi_vals(:);

log_ncgam = -lambda/2 + m + log(sum_a);
end
