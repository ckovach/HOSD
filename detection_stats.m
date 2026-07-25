function [outs,pprobs,xdets,xsnrs,xrsms,xfsds,xfilts,xthrs,gs,pprobs_outlier] = detection_stats(hos, x, varargin)

%[outs,pprob,xdet,xsnr,xrsm] = detection_stats(hos, x, 'noise_dist','gamma', ...)
% Estimates a posterior probability on events detected by xdetect(hos,x)
% under a 2-component mixture model.
%
% INPUTS (required positional):
%   hos - fitted hosobject
%   x   - input signal
%
% INPUTS (optional, positional+name-value):
%   pow - Lp norm exponent for the smoothed-MS statistic (default 2 -> RMS).
%         Accepted as a 3rd positional arg for backward compatibility, or
%         as 'pow', val.
%
% INPUTS (name-value pairs):
%   'noise_dist'   {'gamma','lognormal','weibull','chi2'} (default 'lognormal')
%       Marginal noise family on the smoothed xrsm:
%         'lognormal' : X = xrsm  ~ LogN(mu, sigma). Heavier right tail.
%                       Default based on superiority of fit in the MASS cohort  
%         'gamma'     : Z = xrsm^2 ~ Gamma(a, b). Previous default.
%         'weibull'   : X = xrsm  ~ Weibull(a, b).
%         'chi2'      : Z = xrsm^2 ~ ChiSquared(nu). Single parameter; the
%                       theoretically-correct family under iid Gaussian
%                       xfilt/xfsd. Narrowest of the four.
%   'signal_dist'  {'empirical','gamma','noncentral_gamma','lognormal'} (default 'lognormal')
%       Family of the signal component:
%         'lognormal' : LogN(mu_s, sigma_s). New default
%                       paired with noise_dist='lognormal' + freeze_noise=true.                       
%         'gamma'     : Gamma(a_s, b_s) closed-form M-step. Previous default;
%                       wins F1 vs MODA-exp marginally over the legacy
%                       empirical signal (0.610 vs 0.599) but loses ~0.16
%                       AUPRC vs lognormal on the n=100 cohort.
%			 With freeze_noise=false the joint EM
%                       requires noise_dist='gamma' (gamma noise +
%                       lognormal signal) or noise_dist='lognormal'
%                       (both components lognormal -- routed to
%                       TWO_LOGNORMAL_MIXTURE_EM, with a mu_n < mu_s
%                       label-ordering safety net).
%         'noncentral_gamma' : NCgamma(a_s, b_s, lambda) via Poisson
%                       augmentation; relevant for odd-order HOSD where
%                       the signal need not be zero-mean. Pair with
%                       shared_scale=true to constrain b_n = b_s and let
%                       only lambda differentiate signal from noise. On
%                       symmetric (4th-order) HOSD lambda collapses to 0.
%         'empirical' : Non-parametric histogram on the EM bin grid
%                       (legacy pre-2026-05-20 default; requires
%                       freeze_noise=true). Same F1 as 'gamma' default
%                       (~0.6) but uses 143 free params vs 5; kept for
%                       backwards compatibility / cache regeneration.
%   'freeze_noise' (logical, default true)
%       If true, noise parameters are fit ONCE before the EM (weighted
%       MLE on the kernel-derived sample_w) and held fixed. If false,
%       noise parameters are updated each EM step via responsibilities
%       (joint EM, e.g. the legacy 'two_gamma' string). Default true is
%       conservative -- pinning the noise via sample_w prevents a
%       flexible signal family from cannibalising noise mass.
%   'shared_scale' (logical, default false)
%       Only used when signal_dist='noncentral_gamma' with noise_dist
%       and signal_dist sharing the same family (currently gamma).
%       If true, b_n = b_s = b is enforced -> signal differentiated by
%       lambda alone (4 params total vs 5).
%   'tail_quantile' (default NaN)
%       In (0,1): Peaks-Over-Threshold GPD splice onto the noise pdf
%       (only when freeze_noise=true and signal_dist='empirical'; ignored
%       otherwise). >=50 exceedances required; falls back to bulk-only.
%   'sample_w_pow' (default 3)
%       Exponent on the kernel-derived noise weight sample_w = max(1-wgt,0)^p
%       before the weighted noise MLE. pow=3 lands the gamma fit in a
%       proper unimodal regime (a>1) for all tested MASS subjects and
%       halves the KS-to-noise-CDF vs pow=1. Pass 1 to recover the
%       pre-2026-05-18 behavior.
%   'outlier_trim_quantile' (default NaN)
%       In (0,1): hard quantile of xrsm above which samples are
%       EXCLUDED from the EM fit (their density would otherwise be
%       absorbed by the signal component and bias b_s upward). Samples
%       above the cut are still scored with pprob using the trimmed-fit
%       parameters (they generally come out near 1, as expected for
%       extreme samples under a tight signal model). Currently active
%       only for parametric signal_dist; ignored for empirical.
%
% LEGACY noise_dist values 'two_gamma', 'two_gamma_nc', 'two_gamma_nc_fixed'
% are still accepted (with a one-time deprecation warning) and remapped to:
%   'two_gamma'           -> noise_dist='gamma', signal_dist='gamma',           freeze_noise=false
%   'two_gamma_nc'        -> noise_dist='gamma', signal_dist='noncentral_gamma', freeze_noise=false
%   'two_gamma_nc_fixed'  -> noise_dist='gamma', signal_dist='noncentral_gamma', freeze_noise=false, shared_scale=true
%
% OUTPUTS:
%   outs   - struct array (one element per hos channel) of fit stats. Key
%       distributional parameter fields are grouped into substructs:
%         outs(k).noise_params   .type plus family-specific fitted params
%                                e.g. {gamma:.a,.b,.DF_eff}
%                                     {lognormal:.mu,.sigma}
%                                     {weibull:.a,.b}  {chi2:.nu}
%                                and .mode (x-space mode of fitted density)
%         outs(k).sig_params     analogous, with .type in {gamma,
%                                noncentral_gamma, lognormal, empirical}
%                                noncentral_gamma additionally has .lambda
%         outs(k).outlier_params .type='none' (default) or 'gpd' with
%                                .u, .k, .sigma, .pi, .pi_max
%         outs(k).empirical_pdf  [length(snrx) x 1] empirical density of
%                                the data the EM was fit on, on the same
%                                px0 bin edges as sigpdf/noisepdf. For
%                                visual goodness-of-fit comparison.
%       Per-family parameters that don't apply to the chosen distribution
%       are simply not present in the substruct (no NaN clutter).
%   pprob  - per-sample posterior signal probability
%   xdet,xsnr,xrsm,... - outputs of HOSOBJECT/XDETECT (xrsm/xsnr returned
%       in their THRESHOLDED form; the EM internally uses unthresholded xrsm).
%
%See also HOSOBJECT/XDETECT

%C. Kovach 2026

% ---- Argument parsing ----------------------------------------------------
p = inputParser;
p.addOptional ('pow',                   2,           @(v) isempty(v) || (isnumeric(v) && isscalar(v)));
p.addParameter('noise_dist',            'lognormal', @(s) ischar(s) || isstring(s));
p.addParameter('signal_dist',           'lognormal', @(s) ischar(s) || isstring(s));
p.addParameter('freeze_noise',          [],          @(v) isempty(v) || islogical(v));
p.addParameter('shared_scale',          false,       @islogical);
p.addParameter('tail_quantile',         NaN,         @isnumeric);
p.addParameter('sample_w_pow',          3,           @isnumeric);
p.addParameter('outlier_trim_quantile', NaN,         @isnumeric);
p.addParameter('outlier_dist',          'none',      @(s) ischar(s) || isstring(s));
p.addParameter('outlier_u_quantile',    0.995,       @isnumeric);
p.addParameter('outlier_pi_max',        0.01,        @isnumeric);
% Binning for fast EM: histogram x_fit into K log-spaced bins and run the
% EM on (bin_centers, bin_counts) instead of raw samples. Mathematically
% equivalent to weighted EM up to bin discretization error; at K=500 log
% bins on heavy-tailed z, parameter agreement is ~0.01% with 60-100x
% speedup (see /tmp/verify_binned_em_log.m). Set to 0 to disable.
p.addParameter('bin_em_K',              500,         @isnumeric);
% Force pprob=0 wherever the amplitude statistic xrsm is below the fitted
% noise mode. Below the noise mode the "signal" posterior is only the
% left tail of the signal density overlapping the bulk of the noise, so
% the posterior levitates slightly off zero with no real detection there.
% Flooring it keeps the pprob trace on the floor between events. Default
% true; pass false to recover the raw mixture posterior everywhere.
p.addParameter('floor_below_noise_mode', true,       @islogical);
p.parse(varargin{:});
pow                   = p.Results.pow;        if isempty(pow), pow = 2; end
noise_dist            = lower(char(p.Results.noise_dist));
signal_dist           = lower(char(p.Results.signal_dist));
shared_scale          = p.Results.shared_scale;
tail_quantile         = p.Results.tail_quantile;
sample_w_pow          = p.Results.sample_w_pow;
outlier_trim_quantile = p.Results.outlier_trim_quantile;
outlier_dist          = lower(char(p.Results.outlier_dist));
outlier_u_quantile    = p.Results.outlier_u_quantile;
outlier_pi_max        = p.Results.outlier_pi_max;
floor_below_noise_mode = p.Results.floor_below_noise_mode;
bin_em_K              = p.Results.bin_em_K;
freeze_noise          = p.Results.freeze_noise;

% Legacy two_gamma* string -> map to (noise_dist, signal_dist, freeze_noise, shared_scale)
switch noise_dist
    case 'two_gamma'
        warning('detection_stats:deprecatedTwoGamma', ...
            'noise_dist=''two_gamma'' is deprecated. Use ''signal_dist'',''gamma'',''freeze_noise'',false instead.');
        noise_dist = 'gamma'; signal_dist = 'gamma'; freeze_noise = false;
    case 'two_gamma_nc'
        warning('detection_stats:deprecatedTwoGammaNc', ...
            'noise_dist=''two_gamma_nc'' is deprecated. Use ''signal_dist'',''noncentral_gamma'',''freeze_noise'',false instead.');
        noise_dist = 'gamma'; signal_dist = 'noncentral_gamma'; freeze_noise = false;
    case 'two_gamma_nc_fixed'
        warning('detection_stats:deprecatedTwoGammaNcFixed', ...
            'noise_dist=''two_gamma_nc_fixed'' is deprecated. Use ''signal_dist'',''noncentral_gamma'',''freeze_noise'',false,''shared_scale'',true instead.');
        noise_dist = 'gamma'; signal_dist = 'noncentral_gamma'; freeze_noise = false; shared_scale = true;
end

% Validation
assert(any(strcmp(noise_dist,  {'gamma','lognormal','weibull','chi2'})), ...
    'noise_dist must be one of {gamma, lognormal, weibull, chi2}');
assert(any(strcmp(signal_dist, {'empirical','gamma','noncentral_gamma','lognormal'})), ...
    'signal_dist must be one of {empirical, gamma, noncentral_gamma, lognormal}');
if isempty(freeze_noise)
    freeze_noise = true;   % default: noise pinned via the kernel-weighted MLE.
                           % Pass freeze_noise=false explicitly for the joint
                           % EM (e.g. the legacy two-gamma variants).
end
if strcmp(signal_dist, 'empirical') && ~freeze_noise
    error('detection_stats:badCombination', ...
        ['signal_dist=''empirical'' requires freeze_noise=true.\n', ...
         'An unfrozen empirical signal would absorb arbitrary noise mass (no identifiability).']);
end
if shared_scale && ~strcmp(signal_dist, 'noncentral_gamma')
    warning('detection_stats:sharedScaleIgnored', ...
        'shared_scale only applies to signal_dist=''noncentral_gamma''; ignoring.');
    shared_scale = false;
end
if isfinite(outlier_trim_quantile) && strcmp(signal_dist, 'empirical')
    warning('detection_stats:outlierTrimEmpiricalIgnored', ...
        'outlier_trim_quantile has no effect for signal_dist=''empirical''; ignoring.');
    outlier_trim_quantile = NaN;
end
use_pot   = isfinite(tail_quantile) && tail_quantile > 0 && tail_quantile < 1;
use_trim  = isfinite(outlier_trim_quantile) && outlier_trim_quantile > 0 && outlier_trim_quantile < 1;
parametric_signal = ~strcmp(signal_dist, 'empirical');
assert(any(strcmp(outlier_dist, {'none','gpd'})), ...
    'outlier_dist must be ''none'' or ''gpd''');
use_outlier = strcmp(outlier_dist, 'gpd');
if use_outlier && ~(any(strcmp(signal_dist, {'gamma','lognormal'})) && ~freeze_noise)
    warning('detection_stats:outlierBranchOnly', ...
        ['3-component outlier model is currently implemented only for ', ...
         'signal_dist=''gamma'' or ''lognormal'' with freeze_noise=false. Ignoring outlier_dist.']);
    use_outlier = false;
end

[xdets,xsnrs,xrsms,xfsds,xfilts,xthrs,gs] = xdetect(hos,x,pow);
pprobs_outlier = zeros(size(xrsms));

% %%% Peak weighted xrsms
% xrsms = nthroot(convn((xthrs./xfsds).^pow,gs,'same')./(convn((xthrs>0),gs,'same')+.1),pow);
% xsnrs = xrsms.*xdets;

for k = 1:length(hos)
    xdet = xdets(:,k);
    xsnr = xsnrs(:,k);
    xrsm = xrsms(:,k);
    xfsd = xfsds(:,k);
    xthr = xthrs(:,k);
    xfilt = xfilts(:,k);
    
    g = gs(:,k);

    % --- Empirical noise model: weighted MLE on the smoothed-MS statistic ---
    % xrsm is recomputed from the UNthresholded xfilt (no blanking) so the
    % noise distribution has full empirical support. Sample weights for the
    % noise fit come from the kernel itself:
    %   wgt(t)      = sum_tau g_norm(tau) * 1{xthr(t-tau) ~= 0}
    %               (= kernel-weighted fraction of variance at xrsm(t) that
    %                  comes from supra-threshold/putative-signal samples)
    %   sample_w(t) = max(1 - wgt(t), 0)
    %               (= noise-contribution weight: 1 where the kernel touches
    %                  no event, 0 where the kernel sits entirely on event
    %                  support). Used as MLE weights so every xrsm sample
    %                  contributes to the noise fit in proportion to how
    %                  much of it is genuinely noise-derived. This replaces
    %                  the older hard "dilated mask" (truncated the right
    %                  tail) and "point-only mask" (signal contamination).
    % Smoothed power from the rectified filter output. Handling differs
    % by order parity:
    %   EVEN orders have an unsigned detection statistic (symmetric
    %   threshold [-thr, thr]), so z=|xfilt| is positive almost
    %   everywhere -- plain smoothed power of |xfilt|.
    %   ODD orders have a SIGNED statistic -- the feature has a skewness
    %   direction and the threshold is one-sided ((-Inf, thr], see
    %   @hosobject/xdetect) -- so only positive excursions are genuine
    %   detections. Half-wave via the I(z>0) indicator and average over
    %   the ACTIVE samples only:
    %       xrsm^pow = [g * (I(z>0) .* z^pow)] / [g * I(z>0) + eps]
    %   Discounting the sub-threshold samples (rather than letting them
    %   dilute the average, which full-wave |xfilt| did) keeps a brief
    %   but strong positive transient -- e.g. a K-complex, which is what
    %   the odd/bispectral detector actually picks up -- from being
    %   washed out by the surrounding quiet samples; it is also invariant
    %   to the kernel's overall scale. The active-sample normalization is
    %   skipped for even orders, where the denominator is ~1 anyway.
    if mod(hos(k).order, 2) ~= 0
        z    = xfilt ./ xfsd;   % signed; positive excursions are detections
        pos  = z > 0;
        xrsm = nthroot(convn(pos .* z.^pow, g, 'same') ...
                       ./ (convn(double(pos), g, 'same') + eps), pow);
    else
        xrsm = nthroot(convn(abs(xfilt ./ xfsd).^pow, g, 'same'), pow);
    end
    xsnr = xrsm .* (xdet ~= 0);
    event_mask = double(xthr ~= 0);
    g_norm = g(:) / max(sum(g(:)), eps);
    wgt = convn(event_mask, g_norm, 'same');
    % sample_w_pow > 1 sharpens the weighting: samples sitting half in
    % event support, half in noise (sample_w ~ 0.5) get pushed toward 0,
    % while pure-noise samples (sample_w ~ 1) stay near 1. This
    % counteracts the under-counting of event support that comes from
    % xthr zero-crossings inside spindles -- those samples have falsely
    % high (1 - wgt) and bias the noise-gamma scale upward.
    sample_w = max(1 - wgt, 0) .^ sample_w_pow;

    % =====================================================================
    % Joint-EM parametric mixture (freeze_noise=false, parametric signal).
    % Variants enabled by (signal_dist, shared_scale):
    %   signal_dist='gamma'              : noise ~ Gamma(alpha, beta_n),
    %                                      signal ~ Gamma(alpha, beta_s).
    %                                      4 params (shared shape).
    %   signal_dist='noncentral_gamma'   : noise ~ Gamma(alpha, beta_n),
    %                                      signal ~ NCgamma(alpha, beta_s, lambda).
    %                                      5 params (shared shape, free scales +
    %                                      lambda). With shared_scale=true,
    %                                      beta_n = beta_s = beta: 4 params, only
    %                                      lambda differentiates signal from noise.
    % Bypasses the empirical-sigpdf EM entirely; both component densities
    % are closed-form parametric. Useful for symmetric (4th-order) HOSD,
    % AND for odd-order (e.g. bispectral) HOSD where the signal need not
    % have zero mean -> non-central gamma admits the lambda parameter.
    % =====================================================================
    if parametric_signal && ~freeze_noise
        is_nc       = strcmp(signal_dist, 'noncentral_gamma');
        is_logn     = strcmp(signal_dist, 'lognormal');
        % Two-lognormal joint EM: both noise and signal updated via
        % closed-form weighted MLE on log x at each iteration (see
        % TWO_LOGNORMAL_MIXTURE_EM). This is the only currently-
        % supported non-gamma noise family in the joint-EM branch.
        is_two_logn = strcmp(noise_dist, 'lognormal') && is_logn;
        if ~strcmp(noise_dist, 'gamma') && ~is_two_logn
            error('detection_stats:noiseNotGamma', ...
                ['Joint-EM parametric mixture currently requires noise_dist=''gamma'' (got ''%s''), ', ...
                 'except for the special case (noise_dist=''lognormal'', signal_dist=''lognormal'') ', ...
                 'which routes to TWO_LOGNORMAL_MIXTURE_EM. Mixed families across gamma and lognormal ', ...
                 'would need a separate EM implementation.'], noise_dist);
        end
        nc_mode = '';
        if is_nc && shared_scale,    nc_mode = 'fixed_scale'; end
        if is_nc && ~shared_scale,   nc_mode = 'free_scale';  end

        % Exclude the iterz-induced spike at xrsm ~ 0 (artifact, not a
        % real feature -- biases shape parameter < 1 if included).
        XRSM_MIN = 0.1;
        z_all = xrsm.^2;
        keep_fit = xrsm >= XRSM_MIN & isfinite(z_all);
        % Outlier venting: drop extreme samples from the EM fit so they
        % don't bias beta_s (and, via the shared shape, alpha). Trimmed
        % samples are still scored with the resulting parameters; their
        % pprob naturally lands near 1 under a tight signal model.
        trim_thresh = Inf;
        if use_trim
            trim_thresh = quantile(xrsm(keep_fit), outlier_trim_quantile);
            keep_fit = keep_fit & xrsm <= trim_thresh;
        end
        z_fit = z_all(keep_fit);

        % --- Optional GPD outlier component --------------------------------
        % When use_outlier, pre-fit a Generalized Pareto Distribution to the
        % exceedances of xrsm above an empirical quantile u (in x-space).
        % Convert to a z = x^2 log-density via the Jacobian:
        %   f_o(z) = f_o(x) / (2x)  =>  log f_o(z) = log f_o(x) - log(2 sqrt(z))
        % so the GPD enters the z-space EM with the proper change of variable.
        % The GPD parameters are frozen during EM; only pi_o is estimated,
        % capped at outlier_pi_max so a real-signal mode can't be cannibalised.
        outlier_log_pdf_fit   = [];     % z-space (for gamma EM)
        outlier_log_pdf_fit_x = [];     % x-space (for lognormal EM)
        gpd_u = NaN; gpd_k = NaN; gpd_sigma = NaN;
        if use_outlier
            x_for_gpd = xrsm(keep_fit);
            gpd_u = quantile(x_for_gpd, outlier_u_quantile);
            tail_exc = x_for_gpd(x_for_gpd > gpd_u) - gpd_u;
            if numel(tail_exc) >= 50
                try
                    gp_params = gpfit(tail_exc);
                    gpd_k     = gp_params(1);
                    gpd_sigma = gp_params(2);
                    x_fit_vec = sqrt(z_fit);
                    log_f_o_x = -inf(size(x_fit_vec));
                    above = x_fit_vec > gpd_u;
                    log_f_o_x(above) = log(max(gppdf(x_fit_vec(above) - gpd_u, gpd_k, gpd_sigma), eps));
                    outlier_log_pdf_fit_x = log_f_o_x;
                    % Jacobian transform to z-space: log f_o(z) = log f_o(x) - log(2x)
                    outlier_log_pdf_fit   = log_f_o_x - log(max(2*x_fit_vec, eps));
                catch
                    warning('detection_stats:gpdFitFailed', ...
                        'gpfit failed on tail exceedances; falling back to 2-component model.');
                    use_outlier = false;
                end
            else
                warning('detection_stats:gpdInsufficientTail', ...
                    '<50 tail exceedances above outlier_u_quantile=%.4f; falling back to 2-component model.', ...
                    outlier_u_quantile);
                use_outlier = false;
            end
        end
        % -------------------------------------------------------------------

        x_fit = sqrt(z_fit);

        % --- Pre-bin z_fit (and x_fit) for fast weighted EM ---
        % See the frozen-noise branch below for the rationale and validation
        % (~0.03%-level parameter agreement with raw at K=500 log bins).
        % For the GPD outlier component, log f_o is re-evaluated at bin
        % centers from frozen (gpd_u, gpd_k, gpd_sigma), not subsampled.
        if bin_em_K > 0 && numel(z_fit) > 10*bin_em_K
            z_lo = max(min(z_fit), 1e-6);
            z_hi = quantile(z_fit, 0.99999);
            edges_z  = exp(linspace(log(z_lo), log(z_hi), bin_em_K+1));
            edges_z  = [0, edges_z];
            counts_z = histcounts(z_fit, edges_z);
            centers_z = sqrt(edges_z(1:end-1) .* edges_z(2:end));
            centers_z(1) = edges_z(2) / 2;
            over_z = z_fit(z_fit > z_hi);
            if ~isempty(over_z)
                counts_z  = [counts_z, numel(over_z)];
                centers_z = [centers_z, mean(over_z)];
            end
            keep_b   = counts_z > 0 & centers_z > 0;
            em_z_fit = centers_z(keep_b)';
            em_w_em  = double(counts_z(keep_b))';
            em_x_fit_jt = sqrt(em_z_fit);
            % Re-evaluate outlier log-pdf at bin centers (x-space and z-space)
            if use_outlier
                log_f_o_x_b = -inf(size(em_x_fit_jt));
                above_b = em_x_fit_jt > gpd_u;
                log_f_o_x_b(above_b) = log(max(gppdf(em_x_fit_jt(above_b) - gpd_u, gpd_k, gpd_sigma), eps));
                outlier_log_pdf_em_x = log_f_o_x_b;
                outlier_log_pdf_em_z = log_f_o_x_b - log(max(2*em_x_fit_jt, eps));
            else
                outlier_log_pdf_em_x = [];
                outlier_log_pdf_em_z = [];
            end
        else
            em_z_fit             = z_fit;
            em_w_em              = ones(size(z_fit));
            em_x_fit_jt          = x_fit;
            outlier_log_pdf_em_x = outlier_log_pdf_fit_x;
            outlier_log_pdf_em_z = outlier_log_pdf_fit;
        end

        if is_two_logn
            % Both noise and signal lognormal-on-x. Closed-form weighted
            % MLE on log x for both components every iteration.
            % Asymmetric quantile init places noise at the 30th and
            % signal at the 70th log-quantile (the EM helper enforces
            % mu_n < mu_s with a swap-and-relabel safety net).
            init_em = struct('mu_n', quantile(log(x_fit), 0.30), ...
                             'sigma_n', max(std(log(x_fit(x_fit <= median(x_fit)))), 0.3), ...
                             'mu_s', quantile(log(x_fit), 0.70), ...
                             'sigma_s', max(std(log(x_fit(x_fit >= median(x_fit)))), 0.3), ...
                             'pi_s', 0.10);
            em_args = {'max_iter', 300, 'init', init_em, 'sample_w', em_w_em};
            if use_outlier
                em_args = [em_args, {'outlier_log_pdf', outlier_log_pdf_em_x, ...
                                     'outlier_pi_max',  outlier_pi_max}];
            end
            em = two_lognormal_mixture_em(em_x_fit_jt, em_args{:});
            % Populate gamma-family fields with NaN/0 so the shared
            % downstream code paths (output struct, AIC/BIC) stay
            % uniform with the gamma-noise branches.
            em.alpha  = NaN; em.beta_n = NaN; em.beta_s = NaN; em.lambda = 0;
        elseif is_logn
            init_em = struct('alpha', 2.0, 'beta_n', median(z_fit)/3, ...
                             'mu_s', quantile(log(x_fit), 0.6), ...
                             'sigma_s', max(std(log(x_fit(x_fit >= median(x_fit)))), 0.3), ...
                             'pi_s', 0.10);
            em_args = {'max_iter', 300, 'init', init_em, 'sample_w', em_w_em};
            if use_outlier
                em_args = [em_args, {'outlier_log_pdf', outlier_log_pdf_em_x, ...
                                     'outlier_pi_max',  outlier_pi_max}];
            end
            em = gamma_lognormal_mixture_em(em_x_fit_jt, em_args{:});
            % These fields are not native to the lognormal EM; populate so
            % downstream output-struct code is uniform across signal types.
            em.beta_s = NaN; em.lambda = 0;
        elseif ~is_nc
            init_em = struct('alpha', 2.0, 'beta_n', median(z_fit)/3, ...
                             'beta_s', median(z_fit)*3, 'pi_s', 0.10);
            em_args = {'max_iter', 300, 'init', init_em, 'sample_w', em_w_em};
            if use_outlier
                em_args = [em_args, {'outlier_log_pdf', outlier_log_pdf_em_z, ...
                                     'outlier_pi_max',  outlier_pi_max}];
            end
            em = two_gamma_mixture_em(em_z_fit, em_args{:});
            em.lambda = 0;  % central by construction
            em.mu_s = NaN; em.sigma_s = NaN;
        else
            init_em = struct('alpha',  2.0, ...
                             'beta_n', median(z_fit)/3, ...
                             'beta_s', median(z_fit)*3, ...
                             'lambda', 5.0, ...
                             'pi_s',   0.10);
            if strcmp(nc_mode, 'fixed_scale'); init_em.beta_s = init_em.beta_n; end
            em = gamma_noncentral_mixture_em(em_z_fit, ...
                    'mode',     nc_mode, ...
                    'max_iter', 300, ...
                    'init',     init_em, ...
                    'sample_w', em_w_em);
            em.pi_o = 0;
            em.mu_s = NaN; em.sigma_s = NaN;
        end
        if ~isfield(em, 'pi_o'); em.pi_o = 0; end

        psig       = em.pi_s;
        if is_two_logn
            a_noise     = NaN;
            b_noise     = NaN;
            noise_mode  = exp(em.mu_n - em.sigma_n^2);
            a_signal    = NaN;
            b_signal    = NaN;
            lam_signal  = 0;
            signal_mode = exp(em.mu_s - em.sigma_s^2);
        else
            a_noise    = em.alpha;
            b_noise    = em.beta_n;
            noise_mode = sqrt(max(em.beta_n * (2*em.alpha - 1) / 2, 0));
            if is_logn
                a_signal    = NaN;
                b_signal    = NaN;
                lam_signal  = 0;
                signal_mode = exp(em.mu_s - em.sigma_s^2);
            else
                a_signal    = em.alpha;
                b_signal    = em.beta_s;
                lam_signal  = em.lambda;
                % x-space mode formula (exact for central gamma; approximate for NC)
                signal_mode = sqrt(max(em.beta_s * (2*em.alpha + em.lambda - 1) / 2, 0));
            end
        end

        % Per-sample posterior on ALL samples (responsibility for signal).
        % When an outlier component is active, the 3-way log-sum-exp gives
        % the proper posterior; pprob still reports P(signal | x), with the
        % outlier responsibility tracked separately via pprob_outlier.
        z_safe = max(z_all, eps);
        logz   = log(z_safe);
        % Gamma noise log-pdf on z; signal log-pdf depends on signal_dist.
        % All log-pdfs are expressed on z = x^2 so the 3-way mixture is
        % evaluated in z-space. For the lognormal signal, the native pdf
        % is on x; we apply the Jacobian (log_f_s_z = log_f_s_x - log(2x))
        % so the comparison with f_n(z) is correct.
        % Noise log-pdf in z-space. Gamma noise has a closed-form z-space
        % expression; lognormal noise is native in x-space and needs the
        % log(2x) Jacobian to be commensurate with the gamma branches.
        if is_two_logn
            x_safe = sqrt(z_safe);
            log_x_safe = log(max(x_safe, eps));
            log_f_n_x = -log_x_safe - log(em.sigma_n) - 0.5*log(2*pi) ...
                       - (log_x_safe - em.mu_n).^2 / (2*em.sigma_n^2);
            log_f_n   = log_f_n_x - log(max(2*x_safe, eps));
        else
            log_f_n = (em.alpha-1)*logz - z_safe/em.beta_n - em.alpha*log(em.beta_n) - gammaln(em.alpha);
        end
        if is_logn   % covers both is_two_logn and gamma-noise + lognormal-signal
            x_safe = sqrt(z_safe);
            log_f_s_x = -log(max(x_safe, eps)) - log(em.sigma_s) - 0.5*log(2*pi) ...
                       - (log(max(x_safe, eps)) - em.mu_s).^2 / (2*em.sigma_s^2);
            log_f_s   = log_f_s_x - log(max(2*x_safe, eps));
        elseif ~is_nc
            log_f_s = (em.alpha-1)*logz - z_safe/em.beta_s - em.alpha*log(em.beta_s) - gammaln(em.alpha);
        else
            log_f_s = log_ncgamma_pdf(z_safe, em.alpha, em.beta_s, em.lambda);
        end
        log_pi_n = log(max(1 - psig - em.pi_o, eps));
        log_pi_s = log(max(psig,  eps));
        if use_outlier
            x_all_vec = xrsm;
            log_f_o_x_all = -inf(size(x_all_vec));
            above = x_all_vec > gpd_u & isfinite(x_all_vec);
            log_f_o_x_all(above) = log(max(gppdf(x_all_vec(above) - gpd_u, gpd_k, gpd_sigma), eps));
            log_f_o_all = log_f_o_x_all - log(max(2*sqrt(z_safe), eps));
            log_pi_o = log(max(em.pi_o, eps));
            m_log = max(max(log_pi_n + log_f_n, log_pi_s + log_f_s), log_pi_o + log_f_o_all);
            log_denom = m_log + log(exp(log_pi_n + log_f_n - m_log) ...
                                  + exp(log_pi_s + log_f_s - m_log) ...
                                  + exp(log_pi_o + log_f_o_all - m_log));
            pprob_outlier_all = exp(log_pi_o + log_f_o_all - log_denom);
        else
            m_log = max(log_pi_n + log_f_n, log_pi_s + log_f_s);
            log_denom = m_log + log(exp(log_pi_n + log_f_n - m_log) + exp(log_pi_s + log_f_s - m_log));
            pprob_outlier_all = zeros(size(xrsm));
        end
        pprob = exp(log_pi_s + log_f_s - log_denom);

        % Build px grid and discretize the two component pdfs parametrically.
        snrs_all = xrsm(xdet);
        lb = 0; ub = max([snrs_all; 6]);
        px0 = linspace(lb, ub, (ub-lb)*10);
        px  = (px0(1:end-1) + px0(2:end)) / 2;
        dpx = diff(px(1:2));
        % x-space pdfs (Jacobian 2x: z = x^2 for gamma-family components;
        % lognormal is native on x).
        if is_two_logn
            f_X_noise = @(xx) lognpdf(max(xx, eps), em.mu_n, em.sigma_n);
        else
            f_X_noise = @(xx) 2 .* xx .* gampdf(xx.^2, em.alpha, em.beta_n);
        end
        if is_logn   % covers is_two_logn too
            f_X_signal = @(xx) lognpdf(max(xx, eps), em.mu_s, em.sigma_s);
        elseif ~is_nc
            f_X_signal = @(xx) 2 .* xx .* gampdf(xx.^2, em.alpha, em.beta_s);
        else
            f_X_signal = @(xx) 2 .* xx .* exp(log_ncgamma_pdf(max(xx.^2, eps), em.alpha, em.beta_s, em.lambda));
        end
        noisepdf = f_X_noise(px) * dpx;
        sigpdf   = f_X_signal(px) * dpx;
        noisepdf = noisepdf / max(sum(noisepdf), eps);
        sigpdf   = sigpdf   / max(sum(sigpdf),   eps);
        sigpdf0  = sigpdf;

        snrs = snrs_all(snrs_all >= noise_mode);
        if numel(snrs) < 10; snrs = snrs_all; end

        % --- Output struct (mirrors the standard path's fields) ---
        out.times = find(xdet)./hos(k).sampling_rate;
        out.snrs  = xsnr(xdet);
        out.posterior_prob_estimate = pprob(xdet);
        out.starting_rate_estimate  = sum(xdet)./length(x)*hos(k).sampling_rate;
        out.adj_rate_estimate       = sum(pprob(xdet))./length(x)*hos(k).sampling_rate;
        out.raw_N_estimate      = sum(xdet);
        out.adjusted_N_estimate = sum(pprob(xdet));
        out.est_prior   = psig;
        out.sigpdf      = sigpdf;
        out.sigpdf_raw  = sigpdf0;
        out.noisepdf    = noisepdf;
        % Empirical PDF of the data the EM was fit on (post XRSM_MIN +
        % outlier_trim filters). Same px0 bin edges as sigpdf/noisepdf for
        % direct comparison.
        emp_counts = histcounts(sqrt(z_fit), px0);
        out.empirical_pdf = emp_counts(:) / max(sum(emp_counts), eps);
        % Algorithm flags
        out.freeze_noise = freeze_noise;
        out.shared_scale = shared_scale;
        % Noise / signal / outlier parameter substructs (see make_*_params at end)
        if is_two_logn
            out.noise_params = struct('type',  'lognormal', ...
                                      'mu',    em.mu_n, ...
                                      'sigma', em.sigma_n, ...
                                      'mode',  noise_mode);
        else
            out.noise_params = struct('type', noise_dist, ...
                                      'a',    a_noise, ...
                                      'b',    b_noise, ...
                                      'DF_eff', 2*a_noise, ...
                                      'mode', noise_mode);
        end
        if is_logn   % covers is_two_logn too
            out.sig_params = struct('type',  'lognormal', ...
                                    'mu',    em.mu_s, ...
                                    'sigma', em.sigma_s, ...
                                    'mode',  signal_mode);
        elseif is_nc
            out.sig_params = struct('type',   'noncentral_gamma', ...
                                    'a',      a_signal, ...
                                    'b',      b_signal, ...
                                    'lambda', lam_signal, ...
                                    'mode',   signal_mode);
        else
            out.sig_params = struct('type', 'gamma', ...
                                    'a',    a_signal, ...
                                    'b',    b_signal, ...
                                    'mode', signal_mode);
        end
        if use_outlier
            out.outlier_params = struct('type',  'gpd', ...
                                        'u',     gpd_u, ...
                                        'k',     gpd_k, ...
                                        'sigma', gpd_sigma, ...
                                        'pi',    em.pi_o, ...
                                        'pi_max', outlier_pi_max);
        else
            out.outlier_params = struct('type', 'none');
        end
        % EM diagnostics and preprocessing options.
        % For the two-lognormal mixture there is no closed beta ratio;
        % report the squared median-ratio exp(2 (mu_s - mu_n)) as the
        % natural analog (a "scale^2" comparison between the two log-
        % normal components on the original x-axis).
        if is_two_logn
            out.snr2_ratio = exp(2 * (em.mu_s - em.mu_n));
        else
            out.snr2_ratio = em.beta_s / em.beta_n;
        end
        out.em_iters     = em.n_iter;
        out.em_converged = em.converged;
        out.outlier_trim_quantile = outlier_trim_quantile;
        out.outlier_trim_thresh   = trim_thresh;
        out.pprob_outlier         = pprob_outlier_all(xdet);
        out.pot_active = false;
        out.pot_u = NaN; out.pot_quantile = NaN;
        out.pot_gpd_k = NaN; out.pot_gpd_sigma = NaN; out.pot_p_bulk = NaN;
        out.n_sig_peaks_total = numel(snrs_all);
        out.n_sig_peaks_kept  = numel(snrs);
        out.pow              = pow;
        out.sample_w_pow     = sample_w_pow;
        out.snrx             = px;
        out.thresholds       = px0;
        out.est_recall       = 1 - [0, cumsum(sigpdf)] + eps;
        out.est_specificity  = [0, cumsum(noisepdf)];
        out.est_recall_raw   = 1 - [0, cumsum(sigpdf0)];
        est_precision = psig*(1-[0,cumsum(sigpdf)]) ./ ...
                       (psig*(1-[0,cumsum(sigpdf)]) + (1-psig)*(1-[0,cumsum(noisepdf)]));
        dprecision = diff([0, est_precision]);
        dprecision(dprecision < 0) = 0;
        out.est_precision = cumsum(dprecision);
        out.est_F1     = sqrt(out.est_recall .* out.est_precision);
        out.est_AUROC  = nansum((out.est_recall(1:end-1)+out.est_recall(2:end))/2.*diff(out.est_specificity));
        out.est_AUROC_raw = nansum((out.est_recall_raw(1:end-1)+out.est_recall_raw(2:end))/2.*diff(out.est_specificity));
        out.est_AUPRC  = -nansum((out.est_precision(1:end-1)+out.est_precision(2:end))/2.*diff(out.est_recall));
        [out.max_F1, mxi] = max(out.est_F1);
        out.F1_thresh   = px0(mxi);
        % Cumulant proportion same as the other path
        xfc = cumulant(xfilt/xfsd, hos(k).order, [], false);
        thr = nthroot(hos(k).current_threshold(xfilt/xfsd), hos(k).order);
        if mod(hos(k).order, 2) == 0
            out.cumulant_proportion = 1 + trnormcum(-thr, thr, hos(k).order)./xfc;
        else
            out.cumulant_proportion = 1 + trnormcum(-Inf, thr, hos(k).order)./xfc;
        end
        % GoF: set the frozen-noise-fit-specific diagnostics to NaN; report
        % effective N from the (filtered) data the mixture saw.
        out.gof_ks_noise         = NaN;
        out.gof_tail_95_residual = NaN;
        out.gof_tail_99_residual = NaN;
        out.gof_chi2_noise       = NaN;
        out.gof_chi2_dof         = NaN;
        out.gof_eff_n_noise      = sum(keep_fit);
        out.gof_kl_mixture       = NaN;
        out.gof_chi2_mixture     = NaN;

        % --- AIC / BIC for model comparison (two_gamma family) ---
        % Mixture density evaluated on the valid xrsm samples; when an
        % outlier component is active it adds pi_o * f_o(x) to the density.
        valid_l = isfinite(xrsm) & xrsm > 0;
        nL = sum(valid_l);
        f_n_x = f_X_noise(xrsm(valid_l));
        f_s_x = f_X_signal(xrsm(valid_l));
        if use_outlier
            f_o_x = zeros(nL, 1);
            x_v   = xrsm(valid_l);
            above = x_v > gpd_u;
            f_o_x(above) = gppdf(x_v(above) - gpd_u, gpd_k, gpd_sigma);
            mix_density_l = max((1 - psig - em.pi_o)*f_n_x ...
                              + psig*f_s_x ...
                              + em.pi_o*f_o_x, eps);
        else
            mix_density_l = max(psig*f_s_x + (1-psig)*f_n_x, eps);
        end
        out.logL = sum(log(mix_density_l));
        % Parameter count for AIC/BIC:
        %   gamma signal              : alpha, beta_n, beta_s, pi_s              = 4
        %   noncentral_gamma free     : alpha, beta_n, beta_s, lambda, pi_s       = 5
        %   noncentral_gamma shared   : alpha, beta (shared), lambda, pi_s        = 4
        %   lognormal signal          : alpha, beta_n, mu_s, sigma_s, pi_s        = 5
        % +4 when a GPD outlier component is active (u, k, sigma, pi_o).
        if is_two_logn
            out.n_params = 5;   % mu_n, sigma_n, mu_s, sigma_s, pi_s
        elseif strcmp(signal_dist, 'gamma')
            out.n_params = 4;
        elseif strcmp(signal_dist, 'lognormal')
            out.n_params = 5;
        else  % noncentral_gamma
            if shared_scale, out.n_params = 4; else, out.n_params = 5; end
        end
        if use_outlier
            out.n_params = out.n_params + 4;   % u, k, sigma, pi_o
        end
        out.n_samples_logL = nL;
        out.aic = -2*out.logL + 2*out.n_params;
        out.bic = -2*out.logL + out.n_params*log(nL);

        % Gap fix: where xrsm = 0 (iterz NaN fill-ins propagated through
        % xfilt) both noise and signal pdfs vanish; the log-sum-exp would
        % collapse to the prior pi_s (or pi_s + pi_o under the outlier
        % branch), which downstream consumers read as a real high-posterior
        % event. There is no signal support at a gap -- force pprob = 0.
        gap_mask = ~isfinite(xrsm) | xrsm < 1e-6 ...
                 | (floor_below_noise_mode & xrsm < noise_mode);
        pprob(gap_mask) = 0;
        pprob_outlier_all(gap_mask) = 0;

        outs(k) = out;
        pprobs(:, k) = pprob;
        pprobs_outlier(:, k) = pprob_outlier_all;
        continue;
    end

    % =====================================================================
    % Frozen-noise + parametric-signal branch
    %   noise_dist  in {gamma, lognormal, weibull, chi2}  (any of the 4)
    %   signal_dist in {gamma, noncentral_gamma}
    %   freeze_noise = true
    % Noise parameters come from the weighted MLE on sample_w (same as the
    % empirical-signal frozen path); then a small EM over (signal_params,
    % pi_s) with the noise log-pdf held fixed. Conservative when the
    % signal is rare AND its parametric family is known; AIC/BIC are
    % directly comparable across all three branches because the logL is
    % computed on the same xrsm samples.
    % =====================================================================
    if parametric_signal && freeze_noise
        fmins_opts = optimset('Display','off','TolX',1e-5,'TolFun',1e-6);
        noise_peaks = xrsm;
        keep = noise_peaks > 1e-3 & isfinite(noise_peaks);
        noise_peaks_fit = noise_peaks(keep);
        sample_w_init   = sample_w(keep);
        [a_noise, b_noise, mu_noise, sigma_noise, a_wbl, b_wbl, nu_chi2, ...
            noise_mode, F_X_bulk, f_X_bulk] = ...
            fit_noise_family(noise_peaks_fit, sample_w_init, noise_dist, fmins_opts);

        % Build z (xrsm or xrsm^2 depending on noise family parameterisation)
        % and frozen noise log-pdf at every kept sample.
        XRSM_MIN = 0.1;
        keep_fit = xrsm >= XRSM_MIN & isfinite(xrsm);
        trim_thresh = Inf;
        if use_trim
            trim_thresh = quantile(xrsm(keep_fit), outlier_trim_quantile);
            keep_fit = keep_fit & xrsm <= trim_thresh;
        end
        x_fit  = xrsm(keep_fit);

        % --- Pre-bin x_fit for fast weighted EM ---
        % Histogram into K log-spaced bins; the EM then iterates over ~K
        % weighted "samples" (bin centers, weights = bin counts) instead
        % of the raw ~3M samples. Log spacing resolves the noise mode near
        % zero and the signal mode 10-100x higher; linear bins under-
        % resolve the noise mode and bias the EM (verified against raw EM
        % on 01-02-0001: log-K=200 agrees to 0.02% on all params, raw is
        % the reference). Set bin_em_K=0 to disable and use raw x_fit.
        if bin_em_K > 0 && numel(x_fit) > 10*bin_em_K
            x_lo = max(min(x_fit), 1e-3);
            x_hi = quantile(x_fit, 0.99999);
            edges_x  = exp(linspace(log(x_lo), log(x_hi), bin_em_K+1));
            edges_x  = [0, edges_x];
            counts_x = histcounts(x_fit, edges_x);
            centers_x = sqrt(edges_x(1:end-1) .* edges_x(2:end));
            centers_x(1) = edges_x(2) / 2;  % geometric center for [0, x_lo]
            over_x = x_fit(x_fit > x_hi);
            if ~isempty(over_x)
                counts_x  = [counts_x, numel(over_x)];
                centers_x = [centers_x, mean(over_x)];
            end
            keep_b = counts_x > 0 & centers_x > 0;
            em_x_fit = centers_x(keep_b)';
            em_w     = double(counts_x(keep_b))';
        else
            em_x_fit = x_fit;
            em_w     = ones(size(x_fit));
        end
        W_total = sum(em_w);
        log_f_n_fit = log(max(f_X_bulk(em_x_fit), eps));

        % --- Signal EM over (a_s, b_s, [lambda,] pi_s) with noise frozen.
        % Use the same Poisson-augmentation EM helper as the joint branch,
        % but pass log_f_n as an external 'log_noise_pdf' so the noise
        % parameters are not updated. The helper supports this via the
        % 'frozen_noise_log_pdf' option (see gamma_noncentral_mixture_em.m).
        if strcmp(signal_dist, 'gamma')
            % Closed-form EM: signal-only weighted gamma MLE on responsibilities.
            % We don't have a shared-shape constraint with the noise (noise
            % family may not be gamma), so signal alpha is free.
            a_s = 2.0;
            b_s = max(median(em_x_fit.^2)*2, eps);
            % Same data-driven signal-prior init used in the lognormal
            % branch (see comment there for the audit rationale).
            pi_s = min(max(1 - mean(sample_w), 1e-3), 0.95);
            logL_prev = -Inf;
            for it = 1:200
                % E-step: responsibilities
                log_f_s = log(max(2*em_x_fit .* gampdf(em_x_fit.^2, a_s, b_s), eps));
                log_pi_n = log(max(1-pi_s, eps));
                log_pi_s = log(max(pi_s,   eps));
                m = max(log_pi_n + log_f_n_fit, log_pi_s + log_f_s);
                lden = m + log(exp(log_pi_n + log_f_n_fit - m) + exp(log_pi_s + log_f_s - m));
                gamma_t = exp(log_pi_s + log_f_s - lden);
                logL = sum(em_w .* lden);
                W_s = sum(em_w .* gamma_t);
                pi_s = W_s / W_total;
                % Weighted gamma MLE on z = em_x_fit^2 with weights em_w .* gamma_t
                if W_s > 1
                    z_eff   = em_x_fit.^2;
                    mu_z    = sum(em_w .* gamma_t .* z_eff) / W_s;
                    logmu_z = log(mu_z);
                    mlog_z  = sum(em_w .* gamma_t .* log(max(z_eff, eps))) / W_s;
                    s_ = logmu_z - mlog_z;
                    % Initial via approximation, then 1 Newton step
                    a_s = (3 - s_ + sqrt((s_-3)^2 + 24*s_)) / (12*s_);
                    for ni = 1:50
                        f  = log(a_s) - psi(a_s) - s_;
                        fp = 1/a_s - psi(1, a_s);
                        a_new = a_s - f/fp;
                        if abs(a_new - a_s) < 1e-6*abs(a_s)+1e-9, a_s = a_new; break, end
                        a_s = max(a_new, 1e-3);
                    end
                    b_s = mu_z / a_s;
                end
                if it>1 && abs(logL - logL_prev) < 1e-6*max(abs(logL), 1), break, end
                logL_prev = logL;
            end
            a_signal = a_s; b_signal = b_s; lam_signal = 0;
            mu_s_logn = NaN; sigma_s_logn = NaN;
            em.n_iter = it; em.converged = (it < 200);
            f_X_signal = @(xx) 2 .* xx .* gampdf(xx.^2, a_s, b_s);
            signal_mode = sqrt(max(b_signal * (2*a_signal + lam_signal - 1) / 2, 0));
        elseif strcmp(signal_dist, 'lognormal')
            % Closed-form EM: signal-only weighted lognormal MLE on log x with
            % weights em_w .* gamma_t; the frozen noise log-pdf is held fixed
            % at log_f_n_fit (in x-space). Uses the same log-binned x_fit
            % grid as the gamma-signal branch above.
            logx_fit = log(max(em_x_fit, eps));
            % M-step-first init: gamma_t^(0) = 1 - sample_w is the
            % cumulant-derived signal-responsibility map already used
            % to weight the noise fit. A single weighted lognormal MLE
            % on the raw (un-binned) keep-set with weights gamma_t^(0)
            % supplies (pi_s, mu_s, sigma_s) in one shot, with no
            % heuristic quantile and no lower bound. Selected 2026-06-09
            % after a 100-subject MASS audit (diag_mstep_first.m,
            % diag_mstep_first_moda.m): median +1.37 logL improvement
            % vs the older quantile heuristic, and ZERO subjects with
            % measurable change in AUPRC against MODA gold-standard
            % spindle annotations (median |dAUPRC| < 1e-4). The worst
            % logL "losses" of the new init are visible misfits in
            % the older one where the lognormal-signal component had
            % stretched to cover the noise-mode bulk (see
            % spindle_checker/qm_basins_*.png). Cumulant-consistency
            % parallels the pi_s init audit (reference_detection_stats_pi_init).
            g0_init = max(1 - sample_w(keep_fit), 0);
            log_xkeep = log(max(x_fit, eps));
            W0 = sum(g0_init);
            if W0 > 0
                mu_s    = sum(g0_init .* log_xkeep) / W0;
                sigma_s = sqrt(max(sum(g0_init .* (log_xkeep - mu_s).^2) / W0, eps));
                pi_s    = min(max(mean(g0_init), 1e-3), 0.95);
            else
                cumw     = cumsum(em_w) / W_total;
                mu_s     = logx_fit(find(cumw >= 0.6, 1, 'first'));
                sigma_s  = 0.3;
                pi_s     = min(max(1 - mean(sample_w), 1e-3), 0.95);
            end
            logL_prev = -Inf;
            for it = 1:200
                log_f_s = -logx_fit - log(sigma_s) - 0.5*log(2*pi) - (logx_fit - mu_s).^2 / (2*sigma_s^2);
                log_pi_n = log(max(1-pi_s, eps));
                log_pi_s = log(max(pi_s,   eps));
                m = max(log_pi_n + log_f_n_fit, log_pi_s + log_f_s);
                lden = m + log(exp(log_pi_n + log_f_n_fit - m) + exp(log_pi_s + log_f_s - m));
                gamma_t = exp(log_pi_s + log_f_s - lden);
                logL = sum(em_w .* lden);
                W_s = sum(em_w .* gamma_t);
                pi_s = W_s / W_total;
                if W_s > 1
                    mu_s    = sum(em_w .* gamma_t .* logx_fit) / W_s;
                    sigma_s = sqrt(max(sum(em_w .* gamma_t .* (logx_fit - mu_s).^2) / W_s, 1e-6));
                end
                if it>1 && abs(logL - logL_prev) < 1e-6*max(abs(logL), 1), break, end
                logL_prev = logL;
            end
            a_signal = NaN; b_signal = NaN; lam_signal = 0;
            mu_s_logn = mu_s; sigma_s_logn = sigma_s;
            em.n_iter = it; em.converged = (it < 200);
            f_X_signal = @(xx) lognpdf(max(xx, eps), mu_s, sigma_s);
            signal_mode = exp(mu_s - sigma_s^2);
        else
            error('detection_stats:nyiFrozenSignalNC', ...
                ['signal_dist=''%s'' with freeze_noise=true is not yet implemented. ', ...
                 'Use freeze_noise=false (joint EM) for now.'], signal_dist);
        end

        % Per-sample posterior on ALL xrsm samples using frozen params.
        x_all_safe = max(xrsm, eps);
        log_f_n_all = log(max(f_X_bulk(x_all_safe), eps));
        log_f_s_all = log(max(f_X_signal(x_all_safe), eps));
        log_pi_n = log(max(1-pi_s, eps));
        log_pi_s = log(max(pi_s,   eps));
        m_all  = max(log_pi_n + log_f_n_all, log_pi_s + log_f_s_all);
        ld_all = m_all + log(exp(log_pi_n + log_f_n_all - m_all) + exp(log_pi_s + log_f_s_all - m_all));
        pprob  = exp(log_pi_s + log_f_s_all - ld_all);

        % Build px grid for downstream consumers
        snrs_all = xrsm(xdet);
        lb = 0; ub = max([snrs_all; 6]);
        px0 = linspace(lb, ub, max((ub-lb)*10, 100));
        px  = (px0(1:end-1) + px0(2:end))/2;
        dpx = diff(px(1:2));
        noisepdf = f_X_bulk(px) * dpx; noisepdf = noisepdf / max(sum(noisepdf), eps);
        sigpdf   = f_X_signal(px) * dpx; sigpdf   = sigpdf   / max(sum(sigpdf),   eps);
        sigpdf0  = sigpdf;
        snrs = snrs_all(snrs_all >= noise_mode);
        if numel(snrs) < 10, snrs = snrs_all; end

        % --- Output struct ---
        out.times = find(xdet)./hos(k).sampling_rate;
        out.snrs  = xsnr(xdet);
        out.posterior_prob_estimate = pprob(xdet);
        out.starting_rate_estimate  = sum(xdet)./length(x)*hos(k).sampling_rate;
        out.adj_rate_estimate       = sum(pprob(xdet))./length(x)*hos(k).sampling_rate;
        out.raw_N_estimate      = sum(xdet);
        out.adjusted_N_estimate = sum(pprob(xdet));
        out.est_prior   = pi_s;
        out.sigpdf      = sigpdf;
        out.sigpdf_raw  = sigpdf0;
        out.noisepdf    = noisepdf;
        % Empirical PDF of the data the signal EM was fit on (same px0 bin
        % edges as sigpdf/noisepdf for direct comparison).
        emp_counts = histcounts(x_fit, px0);
        out.empirical_pdf = emp_counts(:) / max(sum(emp_counts), eps);
        out.freeze_noise = freeze_noise;
        out.shared_scale = shared_scale;
        out.noise_params = make_noise_params_local(noise_dist, a_noise, b_noise, ...
            mu_noise, sigma_noise, a_wbl, b_wbl, nu_chi2, noise_mode);
        out.sig_params   = make_sig_params_local(signal_dist, a_signal, b_signal, ...
            lam_signal, mu_s_logn, sigma_s_logn, signal_mode);
        out.outlier_params = struct('type', 'none');
        if ~isnan(b_signal) && ~isnan(b_noise)
            out.snr2_ratio = b_signal / max(b_noise, eps);
        else
            out.snr2_ratio = NaN;
        end
        out.em_iters     = em.n_iter;
        out.em_converged = em.converged;
        out.outlier_trim_quantile = outlier_trim_quantile;
        out.outlier_trim_thresh   = trim_thresh;
        out.pprob_outlier         = zeros(sum(xdet),1);
        out.pot_active = false;
        out.pot_u = NaN; out.pot_quantile = NaN;
        out.pot_gpd_k = NaN; out.pot_gpd_sigma = NaN; out.pot_p_bulk = NaN;
        out.n_sig_peaks_total = numel(snrs_all);
        out.n_sig_peaks_kept  = numel(snrs);
        out.pow              = pow;
        out.sample_w_pow     = sample_w_pow;
        out.snrx             = px;
        out.thresholds       = px0;
        out.est_recall       = 1 - [0, cumsum(sigpdf)] + eps;
        out.est_specificity  = [0, cumsum(noisepdf)];
        out.est_recall_raw   = 1 - [0, cumsum(sigpdf0)];
        est_precision = pi_s*(1-[0,cumsum(sigpdf)]) ./ ...
                       (pi_s*(1-[0,cumsum(sigpdf)]) + (1-pi_s)*(1-[0,cumsum(noisepdf)]));
        dprecision = diff([0, est_precision]); dprecision(dprecision<0)=0;
        out.est_precision = cumsum(dprecision);
        out.est_F1     = sqrt(out.est_recall .* out.est_precision);
        out.est_AUROC  = nansum((out.est_recall(1:end-1)+out.est_recall(2:end))/2.*diff(out.est_specificity));
        out.est_AUROC_raw = nansum((out.est_recall_raw(1:end-1)+out.est_recall_raw(2:end))/2.*diff(out.est_specificity));
        out.est_AUPRC  = -nansum((out.est_precision(1:end-1)+out.est_precision(2:end))/2.*diff(out.est_recall));
        [out.max_F1, mxi] = max(out.est_F1); out.F1_thresh = px0(mxi);

        % Cumulant proportion (same as other paths)
        xfc = cumulant(xfilt/xfsd, hos(k).order, [], false);
        thr = nthroot(hos(k).current_threshold(xfilt/xfsd), hos(k).order);
        if mod(hos(k).order, 2) == 0
            out.cumulant_proportion = 1 + trnormcum(-thr, thr, hos(k).order)./xfc;
        else
            out.cumulant_proportion = 1 + trnormcum(-Inf, thr, hos(k).order)./xfc;
        end

        out.gof_ks_noise         = NaN;
        out.gof_tail_95_residual = NaN;
        out.gof_tail_99_residual = NaN;
        out.gof_chi2_noise       = NaN;
        out.gof_chi2_dof         = NaN;
        out.gof_eff_n_noise      = sum(keep_fit);
        out.gof_kl_mixture       = NaN;
        out.gof_chi2_mixture     = NaN;

        % AIC/BIC over all valid xrsm samples (same denominator as other paths)
        valid_l = isfinite(xrsm) & xrsm > 0;
        nL = sum(valid_l);
        mix_density_l = max(pi_s*f_X_signal(xrsm(valid_l)) + ...
                            (1-pi_s)*f_X_bulk(xrsm(valid_l)), eps);
        out.logL = sum(log(mix_density_l));
        switch noise_dist
            case {'gamma','lognormal','weibull'}, n_noise = 2;
            case 'chi2',                            n_noise = 1;
        end
        switch signal_dist
            case 'gamma',     n_signal_par = 2;
            case 'lognormal', n_signal_par = 2;
            otherwise,        n_signal_par = 3;
        end
        out.n_params = n_noise + n_signal_par + 1;   % noise + signal + pi_s
        out.n_samples_logL = nL;
        out.aic = -2*out.logL + 2*out.n_params;
        out.bic = -2*out.logL + out.n_params*log(nL);

        % Gap fix: see joint-EM block for explanation. Frozen-parametric
        % is the worst-affected branch because both f_X_bulk(x) and
        % f_X_signal(x) carry a 2x Jacobian -> both vanish at x=0 -> both
        % clamp to eps -> pprob collapses to pi_s. Zero out at gaps.
        gap_mask = ~isfinite(xrsm) | xrsm < 1e-6 ...
                 | (floor_below_noise_mode & xrsm < noise_mode);
        pprob(gap_mask) = 0;

        outs(k) = out;
        pprobs(:, k) = pprob;
        continue;
    end

    % =====================================================================
    % Frozen-noise + empirical-signal (legacy / default).
    %   noise_dist  in {gamma, lognormal, weibull, chi2}
    %   signal_dist = 'empirical'
    %   freeze_noise = true
    % =====================================================================
    % noise_peaks kept as an ALL-sample alias for POT splice + diagnostics;
    % the bulk fits below use weighted MLE on the full xrsm with sample_w.
    noise_peaks = xrsm;
    fmins_opts = optimset('Display','off','TolX',1e-5,'TolFun',1e-6);

    keep = noise_peaks > 1e-3 & isfinite(noise_peaks);
    noise_peaks_fit = noise_peaks(keep);
    sample_w_init = sample_w(keep);
    % Sample-level weighted MLE for the initial (and final) noise fit.
    % Run once before the EM; the noise parameters are then FROZEN —
    % only the empirical sigpdf is updated inside the EM.
    [a_noise, b_noise, mu_noise, sigma_noise, a_wbl, b_wbl, nu_chi2, ...
        noise_mode, F_X_bulk, f_X_bulk] = ...
        fit_noise_family(noise_peaks_fit, sample_w_init, noise_dist, fmins_opts);

    % Optional Peaks-Over-Threshold (POT) tail splice. When enabled, we
    % fit a Generalized Pareto Distribution to the exceedances above
    % u = quantile(noise_peaks, tail_quantile), and splice it onto the
    % bulk pdf to produce a valid mixture density:
    %   F_X(x) = (p_bulk / F_bulk(u)) * F_bulk(x)               , x <= u
    %          = p_bulk + (1 - p_bulk) * gpcdf(x - u, k, sigma) , x  > u
    %   f_X(x) = (p_bulk / F_bulk(u)) * f_bulk(x)               , x <= u
    %          = (1 - p_bulk) * gppdf(x - u, k, sigma)          , x  > u
    % where p_bulk = mean(noise_peaks <= u) is the empirical bulk
    % proportion. The bulk family is renormalized to integrate to
    % p_bulk on [0, u]; the GPD carries the (1 - p_bulk) tail mass.
    % Requires >= 50 exceedances; falls back to bulk-only otherwise.
    u_tail = NaN; k_gp = NaN; sigma_gp = NaN; p_bulk = NaN;
    pot_active = false;
    if use_pot
        u_tail = quantile(noise_peaks, tail_quantile);
        tail_exc = noise_peaks(noise_peaks > u_tail) - u_tail;
        if numel(tail_exc) >= 50
            try
                gp_params = gpfit(tail_exc);
                k_gp = gp_params(1);
                sigma_gp = gp_params(2);
                p_bulk = mean(noise_peaks <= u_tail);
                pot_active = true;
            catch
                pot_active = false; % fall back to bulk-only
            end
        end
    end
    % Build sigpdf from xrsm at detection peaks, trimmed below the
    % noise mode of the initial fit.
    snrs_all = xrsm(xdet);
    snrs = snrs_all(snrs_all >= noise_mode);
    if numel(snrs) < 10
        snrs = snrs_all;
    end
    [srt] = sort(snrs);
    lb = 0;
    ub = max([srt;6]);
    px0 = linspace(lb, ub, (ub-lb)*10);
    sigpdf0 = diff(interp1(srt, (1:length(srt))./(length(srt)+1), px0, 'nearest','extrap'));
    px  = (px0(1:end-1) + px0(2:end)) / 2;
    dpx = diff(px(1:2));
    % Normalize sigpdf0 to a proper bin-mass density (the empirical-CDF diff
    % only sums to (N-1)/(N+1) due to nearest-extrap at the endpoints).
    sigpdf0 = sigpdf0 ./ max(nansum(sigpdf0), eps);
    sigpdf  = sigpdf0;

    % Apply optional POT/GPD tail splice and build composite F_X, f_X.
    if pot_active
        F_bulk_u = max(F_X_bulk(u_tail), eps);
        F_X = @(x) (x <= u_tail) .* (p_bulk .* F_X_bulk(x) ./ F_bulk_u) + ...
                   (x >  u_tail) .* (p_bulk + (1 - p_bulk) .* gpcdf(max(x - u_tail, 0), k_gp, sigma_gp));
        f_X = @(x) (x <= u_tail) .* (p_bulk .* f_X_bulk(x) ./ F_bulk_u) + ...
                   (x >  u_tail) .* ((1 - p_bulk) .* gppdf(max(x - u_tail, 0), k_gp, sigma_gp));
    else
        F_X = F_X_bulk; f_X = f_X_bulk;
    end
    % Marginal noise pdf at bin centers, then per-sample normest on the
    % SAME normalization scale.
    %
    % f_X integrates to 1 on [0, infinity), but the px grid only covers
    % [0, ub] (ub set by max signal SNR / floor of 6). Any noise mass past
    % ub -- substantial when POT/GPD with k_gp > 0 is active, or when the
    % noise gamma has a heavy right tail -- is lost. Without renormalization
    % the mixture Bayes ratio pprob_i = psig sigpdf_i / (psig sigpdf_i +
    % (1-psig) noisepdf_i) is biased upward at high xrsm and the EM
    % converges to a fixed point whose psig has absorbed the missing
    % noise mass. We enforce sum(noisepdf) = 1 and rescale normest by the
    % same factor so the per-sample inner-loop pprob uses commensurate
    % signal/noise densities.
    noisepdf = f_X(px) * dpx;
    noisepdf(~isfinite(noisepdf)) = 0;
    noise_norm = sum(noisepdf);
    if abs(noise_norm - 1) > 1e-3
        warning('detection_stats:noisepdfTrunc', ...
                'noisepdf integrates to %.4f over px grid (missing mass = %.4f); renormalizing. Extend ub or fit truncation?', ...
                noise_norm, 1 - noise_norm);
    end
    noisepdf = noisepdf / max(noise_norm, eps);
    normest  = f_X(xrsm) .* dpx / max(noise_norm, eps);
    normest(~isfinite(normest)) = 0;
    normest = max(normest, 1e-6);
    % EM: with noise pdf frozen at the wgt-weighted MLE, only the
    % empirical sigpdf is updated. Inner loop converges psig; outer
    % loop reweights sigpdf via pdfwgt until psig stops moving.
    dpsig0 = inf;
    psig0 = max(mean(xrsm>.5), 1e-3);
    psig = psig0;
    niter = 0;
    while dpsig0 > 1e-6 && niter < 500
        pdfest = interp1(px, sigpdf, xrsm, 'nearest', 'extrap') + eps;
        dpsig = Inf;
        niter2 = 0;
        while dpsig > 1e-6 && niter2 < 100
            pprob = pdfest .* psig ./ (pdfest .* psig + normest .* (1-psig));
            nupprob = mean(pprob);
            dpsig = abs(psig - nupprob);
            psig = nupprob;
            niter2 = niter2 + 1;
        end
        pprob = pdfest .* psig ./ (pdfest .* psig + normest .* (1-psig));
        pdfwgt = sigpdf .* psig ./ (sigpdf .* psig + noisepdf .* (1-psig));
        pdfwgt(isnan(pdfwgt)) = 0;
        sigpdf = sigpdf0 .* pdfwgt;
        sigpdf = sigpdf ./ nansum(sigpdf);
        nupprob0 = mean(pprob);
        dpsig0 = abs(nupprob0 - psig0);
        psig0 = nupprob0;
        psig = psig0;
        niter = niter + 1;
    end
    
    out.times = find(xdet)./hos(k).sampling_rate; %Event times
    out.snrs = xsnr(xdet); %Event snrs
    out.posterior_prob_estimate = pprob(xdet);%Estimated event posterior probability
    out.starting_rate_estimate = sum(xdet)./length(x)*hos(k).sampling_rate;
    out.adj_rate_estimate = sum(pprob(xdet))./length(x)*hos(k).sampling_rate; %Event rates
    
    out.raw_N_estimate= sum(xdet); %Original event number
    out.adjusted_N_estimate = sum(pprob(xdet)); %Estimated event number after accounting for signal-to-noise ratio. If this is in the single digits than the result is doubtful
    out.est_prior = psig; %Estimated prior probability that a sample coincides with an event
    out.sigpdf = sigpdf;
    out.sigpdf_raw = sigpdf0; %Raw distribution of peaksnrs
    out.noisepdf = noisepdf;
    % Empirical PDF on the noise-fit input mask. Same px0 bin edges as
    % sigpdf/noisepdf -- intended for visual goodness-of-fit checking.
    emp_counts = histcounts(noise_peaks_fit, px0);
    out.empirical_pdf = emp_counts(:) / max(sum(emp_counts), eps);
    out.freeze_noise  = freeze_noise;
    out.shared_scale  = shared_scale;
    out.noise_params = make_noise_params_local(noise_dist, a_noise, b_noise, ...
        mu_noise, sigma_noise, a_wbl, b_wbl, nu_chi2, noise_mode);
    % Empirical signal: no parametric params, only the empirical sigpdf above
    out.sig_params   = struct('type', 'empirical');
    out.outlier_params = struct('type', 'none');
    out.outlier_trim_quantile = outlier_trim_quantile;
    out.outlier_trim_thresh   = Inf;       %Only set by parametric-signal branches
    out.pprob_outlier         = zeros(sum(xdet),1);
    out.pot_active    = pot_active;     %true if POT/GPD tail splice was applied
    out.pot_u         = u_tail;         %tail threshold u (xrsm units); NaN if not POT
    out.pot_quantile  = tail_quantile;  %quantile defining u; NaN if not POT
    out.pot_gpd_k     = k_gp;           %GPD shape; NaN if not POT
    out.pot_gpd_sigma = sigma_gp;       %GPD scale; NaN if not POT
    out.pot_p_bulk    = p_bulk;         %empirical bulk proportion below u; NaN if not POT
    out.n_sig_peaks_total = numel(snrs_all);
    out.n_sig_peaks_kept  = numel(snrs);
    out.pow = pow;                 %Lp norm used for xrsm
    out.sample_w_pow = sample_w_pow; %Exponent applied to sample_w
    out.snrx = px; %points at which the pdfs are sampled

    % --- Goodness-of-fit diagnostics on the frozen noise model ---
    % Restrict to the actual MLE input (positive, finite xrsm samples with
    % nontrivial noise weight). The empirical CDF is sample_w-weighted so it
    % targets the same distribution f_X was fit to.
    gof_mask = noise_peaks_fit > 0 & sample_w_init > 1e-6;
    if any(gof_mask)
        x_g = noise_peaks_fit(gof_mask);
        w_g = sample_w_init(gof_mask);
        w_g = w_g(:) / max(sum(w_g(:)), eps);

        % Weighted empirical CDF at the sample points.
        [x_sort, ord] = sort(x_g(:));
        w_sort = w_g(ord);
        F_emp  = cumsum(w_sort);
        F_thy  = F_X(x_sort);

        % (1) Weighted Kolmogorov-Smirnov distance to the fitted noise CDF.
        out.gof_ks_noise = max(abs(F_emp - F_thy));

        % (2) Right-tail residual at the 95th weighted quantile of xrsm.
        %     Theoretical: 1 - F_X(x95) should equal 0.05. Empirical: 0.05.
        %     The discrepancy is 1 - F_X(x95_emp) - 0.05.
        i95 = find(F_emp >= 0.95, 1, 'first');
        if ~isempty(i95)
            x95 = x_sort(i95);
            out.gof_tail_95_residual = (1 - F_X(x95)) - 0.05;
        else
            out.gof_tail_95_residual = NaN;
        end

        % (3) Right-tail residual at the 99th quantile (sensitive to heavy
        %     tails the bulk parametric family underestimates).
        i99 = find(F_emp >= 0.99, 1, 'first');
        if ~isempty(i99)
            x99 = x_sort(i99);
            out.gof_tail_99_residual = (1 - F_X(x99)) - 0.01;
        else
            out.gof_tail_99_residual = NaN;
        end

        % (4) Chi-square GoF on the bulk noise pdf: bin weighted samples on
        %     the EM grid, compare to expected weighted counts under the
        %     fitted noise pdf. Use UN-normalized sample weights so the
        %     absolute total N_w = sum(sample_w) makes the cell-size rule
        %     (exp_n > 5) meaningful. (Earlier version used the
        %     unit-normalized w_g, which made exp_n sub-unity everywhere
        %     and the test reported NaN.)
        sw_g_raw = sample_w_init(gof_mask);
        sw_sort_raw = sw_g_raw(ord);
        N_w = sum(sw_g_raw);
        edges = [0, px0(2:end)];
        [~, bin_idx] = histc(x_sort, edges);                                %#ok<HISTC>
        keep_b = bin_idx >= 1 & bin_idx <= numel(noisepdf);
        if any(keep_b) && N_w > 0
            obs_n = accumarray(bin_idx(keep_b), sw_sort_raw(keep_b), [numel(noisepdf), 1]);
            exp_n = N_w * noisepdf(:);
            ok = exp_n > 5;                          % usual chi2 cell-size rule
            if sum(ok) >= 2
                out.gof_chi2_noise   = sum((obs_n(ok) - exp_n(ok)).^2 ./ exp_n(ok));
                out.gof_chi2_dof     = sum(ok) - 1;  % param count not subtracted (conservative)
            else
                out.gof_chi2_noise = NaN; out.gof_chi2_dof = NaN;
            end
        else
            out.gof_chi2_noise = NaN; out.gof_chi2_dof = NaN;
        end

        % (5) Effective sample size of the weighted noise MLE.
        sw_raw = sample_w_init(:);
        out.gof_eff_n_noise = (sum(sw_raw))^2 / max(sum(sw_raw.^2), eps);
    else
        out.gof_ks_noise         = NaN;
        out.gof_tail_95_residual = NaN;
        out.gof_tail_99_residual = NaN;
        out.gof_chi2_noise       = NaN;
        out.gof_chi2_dof         = NaN;
        out.gof_eff_n_noise      = NaN;
    end

    % (6) Mixture-level GoF: compare empirical histogram of xrsm (ALL
    %     samples, not just noise-weighted) to the EM mixture pdf
    %     psig*sigpdf + (1-psig)*noisepdf on the bin grid. KL divergence
    %     (emp || model) and chi-square.
    [emp_counts, ~] = histcounts(xrsm, [0, px0(2:end)]);
    emp_counts = emp_counts(:);
    n_total    = sum(emp_counts);
    if n_total > 0 && numel(emp_counts) == numel(noisepdf)
        emp_p   = emp_counts / n_total;
        model_p = psig * sigpdf(:) + (1 - psig) * noisepdf(:);
        model_p = model_p / max(sum(model_p), eps);
        nonz    = emp_p > 0 & model_p > 0;
        out.gof_kl_mixture   = sum(emp_p(nonz) .* log(emp_p(nonz) ./ model_p(nonz)));
        out.gof_chi2_mixture = sum((emp_counts(nonz) - n_total*model_p(nonz)).^2 ...
                                   ./ max(n_total*model_p(nonz), eps));
    else
        out.gof_kl_mixture   = NaN;
        out.gof_chi2_mixture = NaN;
    end
    %Predicted performance for various thresholds
    %out.predicted_ppv = sigpdf*rate./(sigpdf*rate+noisepdf*rate);
    out.thresholds = px0;
    out.est_recall = 1-[0,cumsum(sigpdf)]+eps;
    out.est_specificity = [0,cumsum(noisepdf)];
    out.est_recall_raw = 1-[0,cumsum(sigpdf0)];
    est_precision= psig*(1-[0,cumsum(sigpdf)])./(psig*(1-[0,cumsum(sigpdf)])+ (1-psig)*(1-[0,cumsum(noisepdf)]));
    %Force precision to be monotonically increasing, although technically it
    %might not be.
    dprecision = diff([0,est_precision]);
    dprecision(dprecision<0) = 0;
    out.est_precision = cumsum(dprecision);
    out.est_F1 = sqrt(out.est_recall.*out.est_precision);
    out.est_AUROC = nansum((out.est_recall(1:end-1)+out.est_recall(2:end))/2.*diff(out.est_specificity));
    out.est_AUROC_raw = nansum((out.est_recall_raw(1:end-1)+out.est_recall_raw(2:end))/2.*diff(out.est_specificity));
    out.est_AUPRC = -nansum((out.est_precision(1:end-1)+out.est_precision(2:end))/2.*diff(out.est_recall));
    [out.max_F1,mxi] = max(out.est_F1);
    out.F1_thresh = px0(mxi);
    % [out.max_precision,mxi] = max(out.est_precision);
    % out.max_precision_thresh = px0(mxi);
    % probable_N = sum((out.snrs>out.thresholds).*out.est_precision);
    % [out.max_probable_N,mxi] = max(probable_N); 
    
    %%% 1 - of the expected truncated gaussian cumulant to the unthresholded cumulunt
    %%% This is a measure of how much truncation of a Gaussian noise background contributes to
    %%% zeroing out the cumulant. 
    
    xfc = cumulant(xfilt/xfsd,hos(k).order,[],false);
    thr = nthroot(hos(k).current_threshold(xfilt/xfsd),hos(k).order);
    if mod(hos(k).order,2)==0
        out.cumulant_proportion = 1 + trnormcum(-thr,thr,hos(k).order)./xfc;
    else
        out.cumulant_proportion = 1 + trnormcum(-Inf,thr,hos(k).order)./xfc;
    end

    % --- AIC / BIC for model comparison (frozen-noise + empirical-signal) ---
    % logL evaluated as sum of log(mixture density) at each valid xrsm.
    % Signal density = sigpdf(nearest bin)/dpx; noise density = f_X(xrsm)/noise_norm
    % (the renormalisation factor recovers a proper density from the
    % grid-truncated noisepdf). Parameter count is conservative: every
    % bin of the empirical sigpdf is counted as a free parameter (K_bins-1,
    % since it sums to 1), plus the parametric noise's params, pi_s, and
    % +2 if a POT/GPD tail splice is active.
    valid_l = isfinite(xrsm) & xrsm > 0;
    nL = sum(valid_l);
    sigpdf_at = interp1(px, sigpdf, xrsm(valid_l), 'nearest', 'extrap') + eps;
    signal_density_l = sigpdf_at / dpx;
    noise_density_l  = f_X(xrsm(valid_l)) / max(noise_norm, eps);
    mix_density_l = max(psig*signal_density_l + (1-psig)*noise_density_l, eps);
    out.logL = sum(log(mix_density_l));
    switch noise_dist
        case {'gamma','lognormal','weibull'}, n_noise = 2;
        case 'chi2',                            n_noise = 1;
        otherwise,                              n_noise = 2;
    end
    n_signal = numel(px) - 1;     % empirical sigpdf, sums to 1
    n_pot    = 2 * double(pot_active);   % GPD shape, scale
    out.n_params = n_signal + n_noise + 1 + n_pot;   % +1 for pi_s
    out.n_samples_logL = nL;
    out.aic = -2*out.logL + 2*out.n_params;
    out.bic = -2*out.logL + out.n_params*log(nL);

    % Gap fix: see joint-EM block for explanation. The empirical-signal
    % branch already gives pprob ~ 0 at gaps via the normest floor +
    % sigpdf-near-zero combination, but explicitly zeroing here keeps the
    % output semantics uniform across all three branches.
    gap_mask = ~isfinite(xrsm) | xrsm < 1e-6 ...
             | (floor_below_noise_mode & xrsm < noise_mode);
    pprob(gap_mask) = 0;

    outs(k) = out;
    pprobs(:,k) = pprob;
    % Do NOT propagate the unthresholded xrsm/xsnr back into the output
    % arrays. The locally-redefined xrsm (computed from unthresholded xfilt)
    % is what the EM used to produce pprob, but external consumers like
    % get_spindles use the thresholded xrsm from xdetect for morphology
    % peak detection and rely on its peaks matching xdet 1:1. Overwriting
    % the outer xrsms with the unthresholded version produces many extra
    % local maxima and breaks downstream struct-construction.
end
end

% =========================================================================
function [a_n, b_n, mu_n, sig_n, a_w, b_w, nu_c, mode_x, F_X_bulk, f_X_bulk] = ...
    fit_noise_family(x_obs, w_obs, noise_dist, opts)
% Sample-level weighted MLE of the chosen marginal noise family.
% Two of the four families (chi2, lognormal) have closed-form weighted
% MLE; gamma and Weibull use fminsearch. The fit is run ONCE before
% the EM loop and the parameters are frozen — so fminsearch cost is
% paid only once per subject regardless of EM iteration count.
%
% INPUTS
%   x_obs       [N x 1] non-negative xrsm samples (positive)
%   w_obs       [N x 1] non-negative noise-responsibility weights
%                       (typically 1 - wgt, with wgt = kernel-weighted
%                        fraction of variance from supra-threshold xfilt)
%   noise_dist  'gamma' | 'lognormal' | 'weibull' | 'chi2'
%   opts        optimset struct for fminsearch (gamma + Weibull)
%
% OUTPUTS
%   a_n, b_n            gamma shape, scale (NaN unless noise_dist='gamma')
%   mu_n, sig_n         lognormal mu, sigma (NaN unless 'lognormal')
%   a_w, b_w            Weibull scale, shape (NaN unless 'weibull')
%   nu_c                chi2 DF (NaN unless 'chi2')
%   mode_x              mode of the marginal in xrsm units
%   F_X_bulk, f_X_bulk  function handles for the marginal CDF and PDF
a_n = NaN; b_n = NaN; mu_n = NaN; sig_n = NaN;
a_w = NaN; b_w = NaN; nu_c = NaN;
sw = max(sum(w_obs), eps);
% For gamma and Weibull we need numerical MLE; rather than running
% fminsearch on the full sample with continuous weights (slow on
% millions of samples), draw a weighted random subsample and call
% MATLAB's optimized built-in gamfit/wblfit on the equal-weighted
% subsample. Subsampling noise on 100 k samples is well below
% model-misspecification error for our applications.
N_SUB = 100000;
draw_weighted_sub = @(x, w, n) x(randsample(numel(x), min(n, numel(x)), true, max(w, 0) + eps));
switch noise_dist
    case 'gamma'
        try
            x_sub = draw_weighted_sub(x_obs, w_obs, N_SUB);
            gam_params = gamfit(x_sub.^2);
            a_n = gam_params(1);
            b_n = gam_params(2);
        catch
            % method-of-moments fallback
            z = x_obs.^2;
            mu_z  = sum(w_obs .* z) / sw;
            var_z = max(sum(w_obs .* (z - mu_z).^2) / sw, eps);
            a_n = max(mu_z^2 / var_z, 0.1);
            b_n = var_z / mu_z;
        end
        % Mode of f_X(x) = 2x * gampdf(x^2, a, b) in x-space (NOT the bare
        % z-space gamma mode (a-1)*b — the 2x Jacobian shifts it).
        % log f_X = const + (2a-1)log(x) - x^2/b  =>  x_mode = sqrt(b(2a-1)/2)
        % for a > 1/2; mode at 0 otherwise.
        mode_x = sqrt(max(b_n * (2*a_n - 1) / 2, 0));
        F_X_bulk = @(x) gamcdf(x.^2, a_n, b_n);
        f_X_bulk = @(x) 2 .* x .* gampdf(x.^2, a_n, b_n);
    case 'lognormal'
        ly = log(x_obs);
        mu_n  = sum(w_obs .* ly) / sw;
        sig_n = sqrt(max(sum(w_obs .* (ly - mu_n).^2) / sw, 1e-6));
        mode_x = exp(mu_n - sig_n^2);
        F_X_bulk = @(x) logncdf(x, mu_n, sig_n);
        f_X_bulk = @(x) lognpdf(x, mu_n, sig_n);
    case 'chi2'
        z = x_obs.^2;
        nu_c = max(sum(w_obs .* z) / sw, 1);
        % Mode of f_X(x) = 2x * chi2pdf(x^2, nu) — same Jacobian fix as
        % the gamma branch (chi2(nu) = gamma(nu/2, 2), so x_mode = sqrt(nu-1)
        % for nu > 1; mode at 0 otherwise).
        mode_x = sqrt(max(nu_c - 1, 0));
        F_X_bulk = @(x) chi2cdf(x.^2, nu_c);
        f_X_bulk = @(x) 2 .* x .* chi2pdf(x.^2, nu_c);
    case 'weibull'
        try
            x_sub = draw_weighted_sub(x_obs, w_obs, N_SUB);
            wbl_params = wblfit(x_sub);
            a_w = wbl_params(1);
            b_w = wbl_params(2);
        catch
            a_w = max(sum(w_obs .* x_obs) / sw, 1e-3);
            b_w = 2;
        end
        if b_w > 1
            mode_x = a_w * ((b_w - 1) / b_w)^(1/b_w);
        else
            mode_x = 0;
        end
        F_X_bulk = @(x) wblcdf(x, a_w, b_w);
        f_X_bulk = @(x) wblpdf(x, a_w, b_w);
end
end

% =========================================================================
function log_ncgam = log_ncgamma_pdf(z, alpha, beta, lambda, K_max)
% Non-central gamma log-pdf via Poisson-augmentation mixture (stable LSE).
%
%   f_NCgam(z; alpha, beta, lambda) = exp(-lambda/2) * sum_{k>=0}
%      ( (lambda/2)^k / k! ) * Gamma_pdf(z; alpha + k, beta)
%
% INPUTS
%   z      [N x 1] or [1 x N]    non-negative samples
%   alpha  scalar                shape (shared between components)
%   beta   scalar                scale
%   lambda scalar                non-centrality (lambda=0 -> central gamma)
%   K_max  optional scalar       truncation; default
%                                max(5, ceil(lambda/2 + 4*sqrt(lambda/2+1)))
%
% OUTPUT
%   log_ncgam  same shape as z   log f_NCgam evaluated at each z(i)

if nargin < 5 || isempty(K_max)
    K_max = max(5, ceil(lambda/2 + 4 * sqrt(lambda/2 + 1)));
end
sz = size(z);
z  = z(:);
log_z = log(max(z, eps));
log_b = log(max(beta, eps));
log_lh = log(max(lambda/2, eps));
K_vals = 0:K_max;

% log_a_k(z) = k log(lambda/2) - log k! + (alpha+k-1) log z - z/beta
%              - (alpha+k) log beta - log Gamma(alpha+k)
k_only = K_vals * log_lh - gammaln(K_vals + 1) - gammaln(alpha + K_vals);
bracket = K_vals .* (log_z - log_b) + k_only;       % N x (K+1)
common  = (alpha - 1) * log_z - z/beta - alpha*log_b;
log_a   = bracket + common;
m       = max(log_a, [], 2);
sum_a   = sum(exp(log_a - m), 2);
log_ncgam = -lambda/2 + m + log(sum_a);
log_ncgam = reshape(log_ncgam, sz);
end

% =========================================================================
function s = make_noise_params_local(noise_dist, a_n, b_n, mu_n, sigma_n, a_wbl, b_wbl, nu_chi2, noise_mode)
%MAKE_NOISE_PARAMS_LOCAL  Build the out.noise_params substruct.
% Returns a struct with .type plus only the parameters relevant to the
% chosen family, plus .mode. Replaces the older per-family fields
% (out.noise_gamma_a, out.noise_logn_mu, etc.).
s = struct('type', noise_dist);
switch noise_dist
    case 'gamma'
        s.a = a_n; s.b = b_n; s.DF_eff = 2*a_n;
    case 'lognormal'
        s.mu = mu_n; s.sigma = sigma_n;
    case 'weibull'
        s.a = a_wbl; s.b = b_wbl;
    case 'chi2'
        s.nu = nu_chi2;
end
s.mode = noise_mode;
end

function s = make_sig_params_local(signal_dist, a_s, b_s, lambda, mu_s, sigma_s, signal_mode)
%MAKE_SIG_PARAMS_LOCAL  Build the out.sig_params substruct.
s = struct('type', signal_dist);
switch signal_dist
    case 'gamma'
        s.a = a_s; s.b = b_s;
    case 'noncentral_gamma'
        s.a = a_s; s.b = b_s; s.lambda = lambda;
    case 'lognormal'
        s.mu = mu_s; s.sigma = sigma_s;
    case 'empirical'
        % no parametric params
end
s.mode = signal_mode;
end
