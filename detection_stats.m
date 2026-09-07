function [outs,pprobs,xdets,xsnrs,xrsms,xfsds,xfilts,xthrs,gs,pprobs_outlier,local_scales] = detection_stats(hos, x, varargin)

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
%   'options'      struct of any of the parameters below (e.g. the
%                  outs(k).options field of a previous run). Missing fields
%                  take the defaults; explicit name-value pairs given
%                  after the struct override it. The struct may also be
%                  passed as the third positional argument:
%                      detection_stats(hos, x, opts)
%   'local_scale_sec'  running robust noise scale for xfilt (seconds).
%                  0 (default) reproduces the previous behaviour: one
%                  global xfsd for the whole record. A positive value
%                  divides xfilt by a time-varying scale (median |xfilt| in
%                  local_scale_block_sec blocks, supra-threshold samples
%                  excluded, moving median over local_scale_sec, x 1.4826)
%                  before forming xrsm, which removes the scale-mixture
%                  heaviness that non-stationary background produces.
%                  'auto' chooses the span by leave-one-block-out
%                  cross-validation on the block medians: each block's
%                  noise level is predicted from the moving median of the
%                  OTHER blocks in its window and the span with the
%                  smallest squared log prediction error wins (candidates 2..240 s, at least
%                  20 kernel durations and at most a quarter of the
%                  record, plus 'global' = no local scaling). The
%                  candidate set can be given as local_scale_candidates.
%                  Too short a span chases the block noise, too long a span
%                  misses the drift; the CV error is minimal in between,
%                  and if the global scale wins local scaling is switched
%                  off (sec = 0). The scale is returned as local_scales(:,k)
%                  and summarised in outs(k).local_scale (.candidates,
%                  .cv_score show the curve).
%   'local_scale_source'  what the running scale is measured on:
%                  'residual' (default): the component's detection filter
%                  applied to the residual x - sum_j xrec_j, i.e. after
%                  EVERY component's reconstructed signal (thresholded
%                  filter output convolved with the feature waveform,
%                  hosobject/xrec) has been removed. Non-stationarity of
%                  the components themselves - a stimulus-locked event
%                  rate, say - then cannot leak into the noise scale; only
%                  the background does. 'xfilt': the component's own filter
%                  output (which for component k already has components
%                  1..k-1 deflated), with supra-threshold samples excluded.
%                  Supra-threshold samples are excluded from the block
%                  medians in both cases (a safeguard against imperfect
%                  reconstruction).
%   'local_scale_block_sec'  block length of the running scale (default 1)
%   'local_scale_candidates' spans (s) tried by 'auto' (default [] = the
%                  built-in grid [2 3 5 7 10 15 20 30 60 120 240] subject to
%                  the kernel/record limits; a user grid is used as given).
%                  A span shorter than a few stimulus periods will track
%                  stimulus-locked amplitude changes as if they were
%                  background - see outs(k).local_scale.block_acf.
%   'model_select' run several noise/signal/freeze (and optionally
%                  local_scale_sec) combinations and keep, per component,
%                  the best-fitting one. false (default) = off. true =
%                  default grid {gamma,lognormal,weibull} x
%                  {gamma,lognormal}, noise frozen, at the given
%                  local_scale_sec. Or a struct with any of the fields
%                  noise_dist (cellstr), signal_dist (cellstr),
%                  freeze_noise (logical vector), local_scale_sec (numeric
%                  vector); the full factorial grid of valid combinations is
%                  run. Joint-EM candidates (freeze_noise false) are opt-in:
%                  their free 'signal' component tends to absorb the
%                  heteroskedastic bulk of the noise and wins the
%                  likelihood with a prior of 10-20 % and an adjusted N
%                  several times that of the frozen fits, so they should
%                  only be compared against each other. Candidates that
%                  error, whose EM did not converge, whose signal component
%                  lies left of the noise, or whose signal mass is > 5x what
%                  the detections can account for are excluded (in that
%                  order, as long as something remains). The table of
%                  candidates and scores is returned in
%                  outs(k).model_selection.
%   'select_criterion'  'bic' (default) | 'aic' | 'tail' | 'ks'.
%                  Within one local_scale_sec the candidates are compared
%                  on the mixture log-likelihood of the SAME xrsm samples
%                  (BIC/AIC, as already computed by every branch). Across
%                  different local_scale_sec values the likelihoods are not
%                  comparable (the statistic itself changes), so the span
%                  is chosen first by the family-free cross-validation
%                  score described under local_scale_sec = 'auto', and the
%                  family combination is then chosen within that span.
%                  'tail' ranks by |log(observed/nominal exceedance of the
%                  MIXTURE cdf at its 99.9 percent quantile)| and 'ks' by
%                  the weighted KS distance of the noise family (frozen
%                  fits only).
%   'sanity_checks'  'warn' (default) | 'off' | 'error'. Runs the checks
%                  listed under outs(k).checks below and warns (or errors)
%                  when any fails.
%   'noise_dist'   {'gamma','lognormal','weibull','chi2','gengamma'} (default 'lognormal')
%       Marginal noise family on the smoothed xrsm:
%         'lognormal' : X = xrsm  ~ LogN(mu, sigma). Heavier right tail.
%                       Default based on superiority of fit in the MASS cohort  
%         'gamma'     : Z = xrsm^2 ~ Gamma(a, b). Previous default.
%         'weibull'   : X = xrsm  ~ Weibull(a, b).
%         'chi2'      : scaled chi-square, Z = xrsm^2 ~ s * ChiSquared(nu),
%                       nu = 2 E[Z]^2 / Var[Z], s = E[Z]/nu (Satterthwaite
%                       form of an average of n_eff half-normal squares; a
%                       gamma with the shape tied to the variance). Before
%                       2026-09 this family was the unscaled chi2(nu) with
%                       nu = E[Z], which has variance 2 where the normalised
%                       statistic has ~2/n_eff and could not fit.
%         'gengamma'  : X = xrsm ~ generalized gamma (Stacy) with scale a,
%                       shape d and power p: f ~ x^(d-1) exp(-(x/a)^p).
%                       Three parameters; nests gamma-on-x (p=1), Weibull
%                       (d=p), the 'gamma' family above (p=2, d=2a) and the
%                       lognormal as a limit. Weighted MLE via
%                       weighted_gengamma_mle. Frozen fits only (like
%                       weibull/chi2). OPT-IN: not in the default
%                       model_select grid -- pass
%                       model_select=struct('noise_dist',{{'gamma','lognormal','weibull','gengamma'}}).
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
%       outs(k).options        every input parameter as EFFECTIVELY used
%                                (an 'auto' span is replaced by the chosen
%                                span, model_select is off, the selected
%                                families are filled in); pass it back as
%                                'options' to reproduce the fit without
%                                repeating the search. .requested keeps
%                                the original local_scale_sec/model_select.
%       outs(k).local_scale    .source ('residual' | 'xfilt'),
%                                .sec (window used, 0 = global), .p05/.p50/
%                                .p95 and .ratio95to5 of the running scale
%       outs(k).gof_noise_tail_ratio  observed/nominal exceedance of the
%                                fitted noise CDF at the 0.99/0.999/0.9999
%                                quantiles on the noise-weighted xrsm
%                                sample (1 = calibrated, >1 = noise tail too
%                                light, <1 = too heavy/conservative). Note
%                                that residual event mass in the "noise"
%                                sample inflates it.
%       outs(k).gof_mixture_tail_ratio  same for the full mixture CDF on
%                                ALL xrsm samples (the quantity the sanity
%                                check uses)
%       outs(k).checks         sanity checks (.ok plus one logical per
%                                check and .messages): pprob in [0,1] and
%                                finite; posterior non-decreasing in xrsm
%                                above the noise mode; signal component to
%                                the right of the noise component; signal
%                                prior in (0, 0.5); adjusted N <= raw N;
%                                strongest detections classified as signal;
%                                mixture tail calibrated within a factor 5;
%                                EM converged; parameters finite
%       outs(k).model_selection (model_select only) candidate table
%   pprob  - per-sample posterior signal probability
%   xdet,xsnr,xrsm,... - outputs of HOSOBJECT/XDETECT (xrsm/xsnr returned
%       in their THRESHOLDED form; the EM internally uses unthresholded xrsm).
%
%   local_scales - [N x K] running noise scale used for each component
%       (all xfsd when local_scale_sec = 0)
%
%See also HOSOBJECT/XDETECT

%C. Kovach 2026

% ---- Argument parsing ----------------------------------------------------
% An options struct (a previous outs(k).options, or any subset of the
% parameters) may be given as the third positional argument or as
% 'options', S. Its fields are applied first; explicit name-value pairs
% that follow override them. Non-parameter fields are ignored.
[varargin, options_struct] = extract_options_struct(varargin);

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
% Running robust noise scale (see header). 0 = global xfsd (legacy).
p.addParameter('local_scale_sec',       0,           @(v) isnumeric(v) || ischar(v) || isstring(v));
p.addParameter('local_scale_block_sec', 1,           @(v) isnumeric(v) && isscalar(v) && v > 0);
p.addParameter('local_scale_candidates', [],         @(v) isnumeric(v));
p.addParameter('local_scale_source',    'residual',  @(s) ischar(s) || isstring(s));
% Model selection over noise/signal/freeze/local-scale combinations.
p.addParameter('model_select',          false,       @(v) islogical(v) || isstruct(v) || isempty(v));
p.addParameter('select_criterion',      'bic',       @(s) ischar(s) || isstring(s));
% Output sanity checks: 'warn' | 'off' | 'error'.
p.addParameter('sanity_checks',         'warn',      @(v) ischar(v) || isstring(v) || islogical(v));
% Internal: xdetect outputs from a previous call (used by model_select to
% avoid re-running xdetect for every candidate).
p.addParameter('xdetect_cache',         [],          @(v) isempty(v) || isstruct(v));
p.KeepUnmatched = false;
parse_args = merge_option_args(p, options_struct, varargin);
p.parse(parse_args{:});
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
local_scale_sec       = p.Results.local_scale_sec;
local_scale_block_sec = p.Results.local_scale_block_sec;
local_scale_candidates = p.Results.local_scale_candidates(:)';
local_scale_source    = lower(char(p.Results.local_scale_source));
if ~any(strcmp(local_scale_source, {'residual', 'xfilt'}))
    error('detection_stats:badLocalScaleSource', 'local_scale_source must be ''residual'' or ''xfilt''.');
end
model_select          = p.Results.model_select;
select_criterion      = lower(char(p.Results.select_criterion));
sanity_checks         = p.Results.sanity_checks;
xdetect_cache         = p.Results.xdetect_cache;
if islogical(sanity_checks)
    if sanity_checks, sanity_checks = 'warn'; else, sanity_checks = 'off'; end
end
sanity_checks = lower(char(sanity_checks));
assert(any(strcmp(sanity_checks, {'warn','off','error'})), ...
    'sanity_checks must be ''warn'', ''off'' or ''error''');
assert(any(strcmp(select_criterion, {'bic','aic','tail','ks'})), ...
    'select_criterion must be one of {bic, aic, tail, ks}');
if ischar(local_scale_sec) || isstring(local_scale_sec)
    assert(strcmpi(char(local_scale_sec), 'auto'), ...
        'local_scale_sec must be a number of seconds (0 = off) or ''auto''');
    local_scale_sec = 'auto';
else
    assert(isscalar(local_scale_sec) && local_scale_sec >= 0, ...
        'local_scale_sec must be a non-negative scalar or ''auto''');
end

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
assert(any(strcmp(noise_dist,  {'gamma','lognormal','weibull','chi2','gengamma'})), ...
    'noise_dist must be one of {gamma, lognormal, weibull, chi2, gengamma}');
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

% Record every parameter as used (after legacy remapping / validation).
options = struct( ...
    'pow',                    pow, ...
    'noise_dist',             noise_dist, ...
    'signal_dist',            signal_dist, ...
    'freeze_noise',           freeze_noise, ...
    'shared_scale',           shared_scale, ...
    'tail_quantile',          tail_quantile, ...
    'sample_w_pow',           sample_w_pow, ...
    'outlier_trim_quantile',  outlier_trim_quantile, ...
    'outlier_dist',           outlier_dist, ...
    'outlier_u_quantile',     outlier_u_quantile, ...
    'outlier_pi_max',         outlier_pi_max, ...
    'bin_em_K',               bin_em_K, ...
    'floor_below_noise_mode', floor_below_noise_mode, ...
    'local_scale_sec',        local_scale_sec, ...
    'local_scale_block_sec',  local_scale_block_sec, ...
    'local_scale_candidates', local_scale_candidates, ...
    'local_scale_source',     local_scale_source, ...
    'model_select',           model_select, ...
    'select_criterion',       select_criterion, ...
    'sanity_checks',          sanity_checks);

if isempty(xdetect_cache)
    [xdets,xsnrs,xrsms,xfsds,xfilts,xthrs,gs] = xdetect(hos,x,pow);
    scale_src = [];
else
    xdets = xdetect_cache.xdets; xsnrs = xdetect_cache.xsnrs; xrsms = xdetect_cache.xrsms;
    xfsds = xdetect_cache.xfsds; xfilts = xdetect_cache.xfilts; xthrs = xdetect_cache.xthrs;
    gs = xdetect_cache.gs;
    scale_src = xdetect_cache.scale_src;
end
% Signal the running scale is measured on (see 'local_scale_source'). Only
% needed when a local scale is in play; the residual costs one xrec plus
% one filter pass per component.
needs_local = ischar(local_scale_sec) || any(local_scale_sec(:) > 0) || ...
              (~isempty(model_select) && ~(islogical(model_select) && ~model_select));
if needs_local && isempty(scale_src)
    if strcmp(local_scale_source, 'residual')
        [scale_src, ok_res] = residual_filter_output(hos, x, xfilts);
        if ~ok_res
            local_scale_source = 'xfilt'; options.local_scale_source = 'xfilt';
        end
    else
        scale_src = xfilts;
    end
end

% =====================================================================
% Model selection: run every candidate combination (sharing this xdetect
% result), score them, and keep the best per component. Each candidate is
% an ordinary detection_stats call with model_select off.
% =====================================================================
if ~isempty(model_select) && ~(islogical(model_select) && ~model_select)
    cache = struct('xdets', xdets, 'xsnrs', xsnrs, 'xrsms', xrsms, 'xfsds', xfsds, ...
                   'xfilts', xfilts, 'xthrs', xthrs, 'gs', gs, 'scale_src', scale_src);
    [outs, pprobs, pprobs_outlier, local_scales] = ...
        run_model_selection(hos, x, options, model_select, select_criterion, cache);
    for k = 1:numel(outs)
        [outs(k).checks, msgs] = run_sanity_checks(outs(k), pprobs(:,k), xdets(:,k));
        report_sanity(msgs, k, sanity_checks);
    end
    return
end

pprobs_outlier = zeros(size(xrsms));
local_scales   = zeros(size(xfilts));
xrsm_used      = zeros(size(xrsms));
sample_w_used  = zeros(size(xrsms));
local_scale_info = repmat(struct('sec', 0, 'block_sec', local_scale_block_sec, ...
    'p05', NaN, 'p50', NaN, 'p95', NaN, 'ratio95to5', NaN, 'candidates', [], 'cv_score', [], ...
    'cv_se', [], 'block_acf', [], 'block_log_level', [], 'kernel_sec', NaN, ...
    'source', local_scale_source), 1, length(hos));
kernel_sec_used = nan(1, length(hos));

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
    % --- Noise scale: global xfsd (legacy) or a running robust scale ---
    % A slowly drifting background makes the marginal of xrsm a scale
    % mixture that no single family fits in the tail; dividing xfilt by
    % a local scale removes that before the statistic is formed.
    fs_k = hos(k).sampling_rate;
    cv_cands = []; cv_score = []; cv_se = [];
    kern_sec = sum(g > 0) / fs_k;
    kernel_sec_used(k) = kern_sec;
    if isempty(scale_src), src_k = xfilt; else, src_k = scale_src(:,k); end
    if ischar(local_scale_sec)          % 'auto': cross-validated span
        [span_sec, cv_cands, cv_score, cv_se] = choose_local_scale_span(src_k, xthr, fs_k, ...
            local_scale_block_sec, kern_sec, local_scale_candidates);
    else
        span_sec = local_scale_sec;
    end
    if span_sec > 0
        scale_t = running_robust_scale(src_k, xthr, fs_k, span_sec, local_scale_block_sec, xfsd);
    else
        scale_t = xfsd * ones(size(xfilt));
    end
    local_scales(:,k) = scale_t;
    local_scale_info(k).sec = span_sec;
    if ~ischar(local_scale_sec) && span_sec <= 0, local_scale_info(k).source = 'none'; end
    local_scale_info(k).candidates = cv_cands;
    local_scale_info(k).cv_score = cv_score;
    local_scale_info(k).cv_se = cv_se;
    local_scale_info(k).kernel_sec = kern_sec;
    % Autocorrelation of the block noise levels (lags 1..10 blocks): tells
    % whether the level fluctuates on a few-second scale (short spans win
    % the CV) or drifts slowly.
    if span_sec > 0 || ischar(local_scale_sec)
        blk_k = block_medians(src_k, xthr, fs_k, local_scale_block_sec);
        local_scale_info(k).block_acf = block_acf(log(blk_k), 10);
        local_scale_info(k).block_log_level = log(blk_k(:))';
    end
    local_scale_info(k).p05 = prctile(scale_t, 5);
    local_scale_info(k).p50 = prctile(scale_t, 50);
    local_scale_info(k).p95 = prctile(scale_t, 95);
    local_scale_info(k).ratio95to5 = local_scale_info(k).p95 / max(local_scale_info(k).p05, eps);

    if mod(hos(k).order, 2) ~= 0
        z    = xfilt ./ scale_t;   % signed; positive excursions are detections
        pos  = z > 0;
        xrsm = nthroot(convn(pos .* z.^pow, g, 'same') ...
                       ./ (convn(double(pos), g, 'same') + eps), pow);
    else
        xrsm = nthroot(convn(abs(xfilt ./ scale_t).^pow, g, 'same'), pow);
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
    xrsm_used(:,k)     = xrsm;       % kept for the post-loop diagnostics
    sample_w_used(:,k) = sample_w;

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
        if any(strcmp(signal_dist, {'gamma', 'lognormal'}))
            % Signal-only EM with the noise density frozen (helper
            % frozen_signal_em; closed-form weighted MLE M-steps).
            %
            % Initialisation is M-step-first: gamma_t^(0) = 1 - sample_w
            % is the cumulant-derived signal-responsibility map already
            % used to weight the noise fit, and one weighted MLE on the raw
            % keep-set with those weights supplies the signal parameters
            % and pi_s in one shot (selected 2026-06-09 after a 100-subject
            % MASS audit for the lognormal branch: median +1.37 logL vs the
            % older quantile heuristic, zero subjects with measurable AUPRC
            % change against MODA; the gamma branch used a fixed
            % a=2, b=2*median(x^2) start that placed the signal component on
            % the shoulder of the noise bulk and, with a 5-10 % starting
            % prior, could settle there - 640-076 ch73 component 3:
            % prior 0.06, adjusted N 0.1, EM not converged).
            %
            % Degeneracy guard: if the converged signal component's median
            % lies below the noise 95 % quantile, or its prior exceeds 0.25,
            % it is modelling the bulk, not events. The EM is then restarted
            % from a tail init (samples above the noise 99.5 % quantile) and
            % that solution is kept if it is non-degenerate. Which init won
            % is recorded in out.em_init.
            g0_init = max(1 - sample_w(keep_fit), 0);
            W0 = sum(g0_init);
            init1 = struct('pi_s', min(max(1 - mean(sample_w), 1e-3), 0.95));
            if strcmp(signal_dist, 'gamma')
                init1.a = 2.0; init1.b = max(median(em_x_fit.^2) * 2, eps);
                if W0 > 10
                    try
                        [init1.a, init1.b] = weighted_gamma_mle(x_fit.^2, g0_init);
                        init1.pi_s = min(max(mean(g0_init), 1e-3), 0.95);
                    catch
                    end
                end
            else
                log_xkeep = log(max(x_fit, eps));
                if W0 > 0
                    init1.mu    = sum(g0_init .* log_xkeep) / W0;
                    init1.sigma = sqrt(max(sum(g0_init .* (log_xkeep - init1.mu).^2) / W0, eps));
                    init1.pi_s  = min(max(mean(g0_init), 1e-3), 0.95);
                else
                    cumw       = cumsum(em_w) / W_total;
                    init1.mu    = log(max(em_x_fit(find(cumw >= 0.6, 1, 'first')), eps));
                    init1.sigma = 0.3;
                end
            end
            [prm, it, em_conv, ~] = frozen_signal_em(signal_dist, em_x_fit, em_w, log_f_n_fit, init1, 200);
            em_init_used = 'mstep_first';
            q95  = noise_quantile(F_X_bulk, 0.95,  max(x_fit));
            q995 = noise_quantile(F_X_bulk, 0.995, max(x_fit));
            if signal_median_of(signal_dist, prm) < q95 || prm.pi_s > 0.25
                sel = x_fit > q995;
                if sum(sel) >= 10
                    init2 = struct('pi_s', min(max(mean(sel), 1e-3), 0.95));
                    if strcmp(signal_dist, 'gamma')
                        try
                            [init2.a, init2.b] = weighted_gamma_mle(x_fit(sel).^2, ones(sum(sel), 1));
                        catch
                            init2.a = 2.0; init2.b = max(mean(x_fit(sel).^2) / 2, eps);
                        end
                    else
                        lxs = log(max(x_fit(sel), eps));
                        init2.mu = mean(lxs); init2.sigma = max(std(lxs), 0.05);
                    end
                    [prm2, it2, em_conv2, ~] = frozen_signal_em(signal_dist, em_x_fit, em_w, log_f_n_fit, init2, 200);
                    if signal_median_of(signal_dist, prm2) >= q95 && prm2.pi_s <= 0.25
                        prm = prm2; it = it2; em_conv = em_conv2; em_init_used = 'tail';
                    end
                end
            end
            pi_s = prm.pi_s;
            em.n_iter = it; em.converged = em_conv; em.init = em_init_used;
            if strcmp(signal_dist, 'gamma')
                a_s = prm.a; b_s = prm.b;
                a_signal = a_s; b_signal = b_s; lam_signal = 0;
                mu_s_logn = NaN; sigma_s_logn = NaN;
                f_X_signal = @(xx) 2 .* xx .* gampdf(xx.^2, a_s, b_s);
            else
                mu_s = prm.mu; sigma_s = prm.sigma;
                a_signal = NaN; b_signal = NaN; lam_signal = 0;
                mu_s_logn = mu_s; sigma_s_logn = sigma_s;
                f_X_signal = @(xx) lognpdf(max(xx, eps), mu_s, sigma_s);
            end
            signal_mode = signal_mode_of(signal_dist, prm);
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
        if isfield(em, 'init'), out.em_init = em.init; else, out.em_init = ''; end
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
            case 'gengamma',                        n_noise = 3;
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
        case 'gengamma',                        n_noise = 3;
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

% ---- Post-loop diagnostics common to every branch ----------------------
for k = 1:numel(outs)
    % Effective parameters: an 'auto' span is replaced by the span that was
    % chosen, so passing outs(k).options back reproduces this fit without
    % repeating the search. The original request is kept in .requested.
    o_opt = options;
    o_opt.requested = struct('local_scale_sec', options.local_scale_sec, 'model_select', options.model_select);
    if ischar(options.local_scale_sec), o_opt.local_scale_sec = local_scale_info(k).sec; end
    outs(k).options     = o_opt;
    outs(k).local_scale = local_scale_info(k);
    outs(k).kernel_sec  = kernel_sec_used(k);
    outs(k).sampling_rate_used = hos(k).sampling_rate;
    outs(k).gof_noise_tail_ratio = noise_tail_ratio(xrsm_used(:,k), sample_w_used(:,k), ...
        outs(k).thresholds, outs(k).noisepdf, [0.99 0.999 0.9999]);
    mixpdf = outs(k).est_prior * outs(k).sigpdf(:) + (1 - outs(k).est_prior) * outs(k).noisepdf(:);
    outs(k).gof_mixture_tail_ratio = noise_tail_ratio(xrsm_used(:,k), ones(size(xrsm_used(:,k))), ...
        outs(k).thresholds, mixpdf, [0.99 0.999 0.9999]);
    outs(k).gof_tail_quantiles = [0.99 0.999 0.9999];
    outs(k).model_selection = [];
    [outs(k).checks, msgs] = run_sanity_checks(outs(k), pprobs(:,k), xdets(:,k));
    report_sanity(msgs, k, sanity_checks);
end
end


% =========================================================================
function [args, S] = extract_options_struct(args)
%EXTRACT_OPTIONS_STRUCT  Pull an options struct out of varargin.
% Accepts detection_stats(hos, x, S, ...) and/or ..., 'options', S, ...
S = struct();
if ~isempty(args) && isstruct(args{1})
    S = args{1}; args = args(2:end);
end
k = 1;
while k <= numel(args) - 1
    if (ischar(args{k}) || isstring(args{k})) && strcmpi(char(args{k}), 'options')
        if isstruct(args{k+1})
            f = fieldnames(args{k+1});
            for j = 1:numel(f), S.(f{j}) = args{k+1}.(f{j}); end
        end
        args(k:k+1) = [];
    else
        k = k + 1;
    end
end
end

function args = merge_option_args(p, S, explicit)
%MERGE_OPTION_ARGS  Struct fields first (known parameters only), then the
% explicit name-value pairs, so that explicit pairs win.
known = p.Parameters;
args = {};
lead = {};
% keep a leading positional pow if one was given explicitly
if ~isempty(explicit) && ~(ischar(explicit{1}) || isstring(explicit{1}))
    lead = explicit(1); explicit = explicit(2:end);
end
f = fieldnames(S);
skipped = {};
for j = 1:numel(f)
    if any(strcmpi(f{j}, known))
        if strcmpi(f{j}, 'pow')
            if isempty(lead), lead = {S.(f{j})}; end
        elseif strcmpi(f{j}, 'xdetect_cache')
            % never carried in an options struct
        else
            args(end+1:end+2) = {f{j}, S.(f{j})};
        end
    else
        skipped{end+1} = f{j}; %#ok<AGROW>
    end
end
if ~isempty(skipped)
    % Output structs carry many non-parameter fields; only complain about
    % names that look like parameters (short, lower-case, no 'est_'/'gof_').
    odd = skipped(cellfun(@(s) isempty(regexp(s, '^(est_|gof_|pot_|n_|em_|requested$)', 'once')) && numel(s) < 24, skipped));
    if ~isempty(odd) && numel(odd) < 8
        warning('detection_stats:unknownOptionFields', ...
            'Ignoring unknown option field(s): %s', strjoin(odd, ', '));
    end
end
args = [lead, args, explicit];
end

function [xsrc, ok] = residual_filter_output(hos, x, xfilts)
%RESIDUAL_FILTER_OUTPUT  Each component's detection filter applied to the
% residual after removing EVERY component's reconstructed signal
% (hosobject/xrec: thresholded filter output * feature waveform, LMSE
% scaled). Falls back to the plain filter outputs if the reconstruction is
% unavailable (e.g. hos is not a hosobject).
xsrc = xfilts; ok = true;
try
    xr = xrec(hos, x(:));                     % N x K (deflated, all components)
    xr(isnan(xr)) = 0;
    xres = x(:) - sum(xr, 2);
    for k = 1:numel(hos)
        xf = apply_filter(hos(k), xres, false, false);
        xsrc(:,k) = xf(:);
    end
catch ME
    ok = false;
    warning('detection_stats:residualScale', ...
        'Could not form the residual filter output (%s); the running scale is measured on xfilt instead.', ME.message);
end
end

function [blk, nb, nBlk] = block_medians(xfilt, xthr, fs, block_sec)
%BLOCK_MEDIANS  Median |xfilt| per block_sec block, supra-threshold samples
% excluded (NaN where a block has no usable sample).
nb   = max(1, round(block_sec * fs));
nBlk = ceil(numel(xfilt) / nb);
xa = abs(xfilt(:)); xa(xthr(:) ~= 0) = NaN;
xa(end+1:nBlk*nb) = NaN;
blk = median(reshape(xa, nb, nBlk), 1, 'omitnan')';
end

function [best, cands, score, se] = choose_local_scale_span(xfilt, xthr, fs, block_sec, kern_sec, user_cands)
%CHOOSE_LOCAL_SCALE_SPAN  Cross-validated span for the running scale.
% Leave-one-block-out CV on the block medians: the log level of each block
% is predicted by the moving median (span) of the OTHER blocks in its
% window; score = mean squared log prediction error, se = its standard
% error. Inf = global scale (median of all other blocks). The candidate
% with the smallest score wins; 0 is returned when the global scale wins
% (local scaling off).
[blk, ~, nBlk] = block_medians(xfilt, xthr, fs, block_sec);
T = nBlk * block_sec;
if nargin < 6 || isempty(user_cands)
    cands = [2 3 5 7 10 15 20 30 60 120 240];
    cands = cands(cands >= 20 * kern_sec & cands <= T / 4);
else
    cands = sort(user_cands(isfinite(user_cands) & user_cands > 0));
end
cands = [cands, Inf];
score = nan(size(cands)); se = nan(size(cands));
for c = 1:numel(cands)
    [score(c), se(c)] = cv_span_score(blk, nBlk, block_sec, cands(c));
end
[~, j] = min(score);
if isempty(j) || ~isfinite(score(j)) || isinf(cands(j))
    best = 0;
else
    best = cands(j);
end
end

function [sc, se] = cv_span_score(blk, nBlk, block_sec, span_sec)
%CV_SPAN_SCORE  Leave-one-block-out squared log prediction error of the
% block levels for one span (0/Inf = global). Returns the mean and its SE.
ok = isfinite(blk) & blk > 0; lb = log(blk); lb(~ok) = NaN;
if span_sec <= 0 || isinf(span_sec)
    h = nBlk;
else
    h = floor(max(3, round(span_sec / block_sec)) / 2);
end
pred = loo_moving_median(lb, h);
e = lb(ok) - pred(ok); e = e(isfinite(e));
if isempty(e), sc = Inf; se = Inf; else, sc = mean(e.^2); se = std(e.^2) / sqrt(numel(e)); end
end

function m = loo_moving_median(v, h)
%LOO_MOVING_MEDIAN  Median of the finite values in window i-h..i+h EXCLUDING
% block i itself (NaN where the window has no other finite value).
v = v(:); N = numel(v); m = nan(N, 1);
for i = 1:N
    lo = max(1, i - h); hi = min(N, i + h);
    w = [v(lo:i-1); v(i+1:hi)]; w = w(isfinite(w));
    if ~isempty(w), m(i) = median(w); end
end
end

function r = block_acf(v, maxlag)
%BLOCK_ACF  Sample autocorrelation of a (NaN-tolerant) block series, lags 1..maxlag.
v = v(:); v = v(isfinite(v)); v = v - mean(v); n = numel(v);
r = nan(1, maxlag);
if n < 3 * maxlag, return, end
den = sum(v.^2);
for L = 1:maxlag
    r(L) = sum(v(1:end-L) .* v(1+L:end)) / max(den, eps);
end
end

function m = moving_median_omitnan(v, n)
%MOVING_MEDIAN_OMITNAN  Centred moving median of length n ignoring NaNs
% (NaN where a window has no finite value). Block-level vectors only.
v = v(:); N = numel(v); h = floor(n / 2); m = nan(N, 1);
for i = 1:N
    w = v(max(1, i - h):min(N, i + h)); w = w(isfinite(w));
    if ~isempty(w), m(i) = median(w); end
end
end

function scale_t = running_robust_scale(xfilt, xthr, fs, span_sec, block_sec, xfsd)
%RUNNING_ROBUST_SCALE  Time-varying noise scale of the filter output.
% Median |xfilt| in block_sec blocks (supra-threshold samples excluded),
% moving median over span_sec, linear interpolation back to samples,
% x 1.4826 (median |x| -> SD for Gaussian noise). Falls back to xfsd
% wherever the estimate is undefined.
[blk, nb, nBlk] = block_medians(xfilt, xthr, fs, block_sec);
nspan = max(3, round(span_sec / block_sec));
blk = moving_median_omitnan(blk, nspan);
bad = ~isfinite(blk) | blk <= 0;
if all(bad)
    blk(:) = xfsd / 1.4826;
elseif any(bad)
    blk(bad) = interp1(find(~bad), blk(~bad), find(bad), 'nearest', 'extrap');
end
tb = ((0:nBlk-1)' + 0.5) * nb;
scale_t = interp1(tb, blk, (1:numel(xfilt))', 'linear', 'extrap') * 1.4826;
scale_t(~isfinite(scale_t) | scale_t <= 0) = xfsd;
scale_t = reshape(scale_t, size(xfilt));
end

function tr = noise_tail_ratio(xrsm, sample_w, px0, noisepdf, qs)
%NOISE_TAIL_RATIO  observed/nominal exceedance of the fitted noise CDF on
% the noise-weighted xrsm sample. noisepdf is the bin-mass vector on the
% px0 edges (as stored in out.noisepdf / out.thresholds).
tr = nan(size(qs));
keep = isfinite(xrsm) & xrsm >= 0.1 & sample_w > 0;
xs = xrsm(keep); ws = sample_w(keep);
if isempty(xs) || sum(ws) <= 0 || isempty(noisepdf), return, end
Fn = [0; cumsum(noisepdf(:))]; Fn = Fn / max(Fn(end), eps);
px0 = px0(:);
if numel(px0) ~= numel(Fn), return, end
for i = 1:numel(qs)
    j = find(Fn >= qs(i), 1, 'first');
    if isempty(j), continue, end
    if j == 1, thr = px0(1);
    else, thr = px0(j-1) + (qs(i) - Fn(j-1)) / max(Fn(j) - Fn(j-1), eps) * (px0(j) - px0(j-1));
    end
    tr(i) = (sum(ws(xs > thr)) / sum(ws)) / (1 - qs(i));
end
end

function [checks, msgs] = run_sanity_checks(out, pprob, xdet)
%RUN_SANITY_CHECKS  Does the fitted mixture make sense?
checks = struct('ok', true);
msgs = {};

% 1. posterior values
bad = ~isfinite(pprob) | pprob < 0 | pprob > 1;
[checks, msgs] = mark(checks, msgs, 'pprob_in_range', ~any(bad), ...
    sprintf('pprob has %d non-finite or out-of-range values', sum(bad)));

% 2. posterior non-decreasing in xrsm above the noise mode (on the model grid)
px = out.snrx(:); ps = out.est_prior;
mix = ps * out.sigpdf(:) + (1 - ps) * out.noisepdf(:);
post = ps * out.sigpdf(:) ./ max(mix, eps);
nm = out.noise_params.mode; if ~isfinite(nm), nm = 0; end
use = find(px >= nm & mix > 1e-4 * max(mix));   % ignore the numerically empty far tail
mono = true; msg2 = '';
if numel(use) >= 3
    d = diff(post(use));
    [dmin, imin] = min(d);
    mono = dmin >= -0.02;
    if ~mono
        msg2 = sprintf(['posterior decreases with xrsm above the noise mode ' ...
            '(largest drop %.3f at xrsm ~ %.2f); check that the signal component lies to the right of the noise'], ...
            -dmin, px(use(imin)));
    end
end
[checks, msgs] = mark(checks, msgs, 'pprob_monotone', mono, msg2);

% 3. signal component to the right of the noise component
mn = sum(px .* out.noisepdf(:)) / max(sum(out.noisepdf), eps);
ms = sum(px .* out.sigpdf(:))   / max(sum(out.sigpdf),   eps);
[checks, msgs] = mark(checks, msgs, 'signal_right_of_noise', ms > mn, ...
    sprintf('signal mean (%.2f) is not above the noise mean (%.2f)', ms, mn));

% 4. signal prior
[checks, msgs] = mark(checks, msgs, 'prior_reasonable', isfinite(ps) && ps > 0 && ps < 0.5, ...
    sprintf('estimated signal prior %.3g is outside (0, 0.5)', ps));

% 5. adjusted N
okN = out.adjusted_N_estimate <= out.raw_N_estimate + 1e-9 && out.adjusted_N_estimate >= 0;
[checks, msgs] = mark(checks, msgs, 'adjustedN_le_rawN', okN, ...
    sprintf('adjusted N (%.1f) exceeds raw N (%d)', out.adjusted_N_estimate, out.raw_N_estimate));

% 6. the strongest detections should be classified as signal
okTop = true; top = NaN;
if any(xdet)
    ppd = pprob(xdet); sn = out.snrs(:);
    [~, ord] = sort(sn, 'descend');
    ntop = max(1, ceil(0.02 * numel(ord)));
    top = mean(ppd(ord(1:ntop)));
    okTop = top >= 0.5;
end
[checks, msgs] = mark(checks, msgs, 'top_detections_are_signal', okTop, ...
    sprintf(['mean posterior of the top 2%% of detections is %.2f: the signal component ' ...
             'does not capture the strongest events (nothing is detected)'], top));

% 7. mixture tail calibration (0.999 quantile of the fitted mixture on all samples)
okTail = true; msg7 = '';
if isfield(out, 'gof_mixture_tail_ratio') && numel(out.gof_mixture_tail_ratio) >= 2
    tr = out.gof_mixture_tail_ratio(2);
    okTail = ~isfinite(tr) || (tr > 1/5 && tr < 5);
    if ~okTail
        trn = NaN; if isfield(out, 'gof_noise_tail_ratio'), trn = out.gof_noise_tail_ratio(2); end
        msg7 = sprintf(['mixture mis-calibrated in the tail: observed/nominal exceedance at its ' ...
            '99.9%% quantile = %.2f (%s; noise-only ratio %.2f)'], tr, ...
            ternary(tr > 1, 'tail too light -> false positives', 'tail too heavy'), trn);
    end
end
[checks, msgs] = mark(checks, msgs, 'mixture_tail_calibrated', okTail, msg7);

% 8. EM convergence
okEM = ~(isfield(out, 'em_converged') && ~isempty(out.em_converged) && ~out.em_converged);
itn = NaN; if isfield(out, 'em_iters'), itn = out.em_iters; end
[checks, msgs] = mark(checks, msgs, 'em_converged', okEM, sprintf('EM did not converge (%d iterations)', itn));

% 9. parameters finite
okp = true;
fn = fieldnames(out.noise_params);
for i = 1:numel(fn)
    v = out.noise_params.(fn{i}); if isnumeric(v) && any(~isfinite(v)), okp = false; end
end
fn = fieldnames(out.sig_params);
for i = 1:numel(fn)
    v = out.sig_params.(fn{i});
    if isnumeric(v) && any(~isfinite(v)) && ~strcmp(fn{i}, 'mode'), okp = false; end
end
[checks, msgs] = mark(checks, msgs, 'params_finite', okp, 'non-finite noise/signal parameters');

% 10. signal-component mass vs what the detections can account for.
% Each detected event raises xrsm over roughly one kernel support, so the
% signal prior should be of the order raw_N * kernel_samples / N. A
% component holding many times that mass is modelling the background
% (typically the joint-EM 'signal' absorbing the heteroskedastic bulk).
okMass = true; ratio = NaN;
if isfield(out, 'kernel_sec') && isfinite(out.kernel_sec) && out.kernel_sec > 0 && out.raw_N_estimate > 0 ...
        && isfield(out, 'sampling_rate_used') && isfinite(out.sampling_rate_used)
    kern_samples = max(1, out.kernel_sec * out.sampling_rate_used);
    ratio = ps * numel(pprob) / (out.raw_N_estimate * kern_samples);
    okMass = ratio <= 5;
end
[checks, msgs] = mark(checks, msgs, 'signal_mass_vs_detections', okMass, ...
    sprintf(['signal component holds %.1f x the mass the %d detections can account for ' ...
             '(prior %.3g); it is probably absorbing background rather than events'], ...
             ratio, out.raw_N_estimate, ps));
checks.messages = msgs;
end

function [checks, msgs] = mark(checks, msgs, name, ok, msg)
checks.(name) = logical(ok);
if ~ok
    checks.ok = false;
    msgs{end+1} = msg;
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function report_sanity(msgs, k, mode)
if isempty(msgs) || strcmp(mode, 'off'), return, end
txt = sprintf('detection_stats sanity checks failed for component %d:\n  - %s', k, strjoin(msgs, sprintf('\n  - ')));
if strcmp(mode, 'error')
    error('detection_stats:sanityCheck', '%s', txt);
else
    warning('detection_stats:sanityCheck', '%s', txt);
end
end

function [outs, pprobs, pprobs_outlier, local_scales] = ...
    run_model_selection(hos, x, options, model_select, criterion, cache)
%RUN_MODEL_SELECTION  Fit every candidate combination and keep the best per
% component.
G = struct('noise_dist', {{'gamma','lognormal','weibull'}}, ...
           'signal_dist', {{'gamma','lognormal'}}, ...
           'freeze_noise', true, ...
           'local_scale_sec', {{options.local_scale_sec}});
if isstruct(model_select)
    if isfield(model_select, 'noise_dist'),   G.noise_dist   = cellstr(model_select.noise_dist);   end
    if isfield(model_select, 'signal_dist'),  G.signal_dist  = cellstr(model_select.signal_dist);  end
    if isfield(model_select, 'freeze_noise'), G.freeze_noise = logical(model_select.freeze_noise(:))'; end
    if isfield(model_select, 'local_scale_sec')
        v = model_select.local_scale_sec;
        if isnumeric(v), G.local_scale_sec = num2cell(v(:))'; else, G.local_scale_sec = cellstr(v); end
    end
end
% Candidate grid (valid combinations only)
C = struct('noise_dist', {}, 'signal_dist', {}, 'freeze_noise', {}, 'local_scale_sec', {});
for ls = 1:numel(G.local_scale_sec)
  for fz = G.freeze_noise
    for nd = G.noise_dist
      for sd = G.signal_dist
        okc = true;
        if strcmp(sd{1}, 'empirical') && ~fz, okc = false; end
        if ~fz && ~(strcmp(nd{1}, 'gamma') || (strcmp(nd{1}, 'lognormal') && strcmp(sd{1}, 'lognormal'))), okc = false; end
        if okc
            C(end+1) = struct('noise_dist', nd{1}, 'signal_dist', sd{1}, 'freeze_noise', fz, ...
                              'local_scale_sec', G.local_scale_sec{ls}); %#ok<AGROW>
        end
      end
    end
  end
end
nC = numel(C); K = numel(hos); N = numel(x);
fprintf('detection_stats model selection: %d candidate(s) x %d component(s)\n', nC, K);
res = cell(1, nC);
for c = 1:nC
    oc = options;
    oc.noise_dist = C(c).noise_dist; oc.signal_dist = C(c).signal_dist;
    oc.freeze_noise = C(c).freeze_noise; oc.local_scale_sec = C(c).local_scale_sec;
    oc.model_select = false; oc.sanity_checks = 'off';
    if ~strcmp(oc.signal_dist, 'noncentral_gamma'), oc.shared_scale = false; end
    try
        [o_c, pp_c, ~, ~, ~, ~, ~, ~, ~, ppo_c, ls_c] = detection_stats(hos, x, oc, 'xdetect_cache', cache);
        res{c} = struct('outs', o_c, 'pprobs', pp_c, 'pprobs_outlier', ppo_c, 'local_scales', ls_c, 'err', '');
    catch ME
        res{c} = struct('outs', [], 'pprobs', [], 'pprobs_outlier', [], 'local_scales', [], 'err', ME.message);
        fprintf('  candidate %d (%s/%s, freeze=%d, scale=%s) failed: %s\n', c, C(c).noise_dist, ...
            C(c).signal_dist, C(c).freeze_noise, num2str(C(c).local_scale_sec), ME.message);
    end
end
% Score table and selection per component
pprobs = zeros(N, K); pprobs_outlier = zeros(N, K); local_scales = zeros(N, K);
outs = [];
for k = 1:K
    T = repmat(struct('noise_dist', '', 'signal_dist', '', 'freeze_noise', false, 'local_scale_sec', 0, ...
        'bic', NaN, 'aic', NaN, 'logL', NaN, 'n_params', NaN, 'tail_ratio_999', NaN, 'noise_tail_ratio_999', NaN, ...
        'ks_noise', NaN, 'cv_span_score', NaN, ...
        'adjusted_N', NaN, 'est_prior', NaN, 'signal_right_of_noise', false, 'em_converged', true, ...
        'mass_ok', true, 'error', '', 'selected', false), 1, nC);
    for c = 1:nC
        T(c).noise_dist = C(c).noise_dist; T(c).signal_dist = C(c).signal_dist;
        T(c).freeze_noise = C(c).freeze_noise; T(c).local_scale_sec = C(c).local_scale_sec;
        T(c).error = res{c}.err;
        if isempty(res{c}.outs), continue, end
        o = res{c}.outs(k);
        T(c).bic = o.bic; T(c).aic = o.aic; T(c).logL = o.logL; T(c).n_params = o.n_params;
        T(c).tail_ratio_999 = o.gof_mixture_tail_ratio(2); T(c).ks_noise = o.gof_ks_noise;
        T(c).noise_tail_ratio_999 = o.gof_noise_tail_ratio(2);
        T(c).cv_span_score = NaN;
        T(c).adjusted_N = o.adjusted_N_estimate; T(c).est_prior = o.est_prior;
        px = o.snrx(:);
        T(c).signal_right_of_noise = sum(px .* o.sigpdf(:)) / max(sum(o.sigpdf), eps) > ...
                                     sum(px .* o.noisepdf(:)) / max(sum(o.noisepdf), eps);
        T(c).em_converged = ~(isfield(o, 'em_converged') && ~isempty(o.em_converged) && ~o.em_converged);
        if isfinite(o.kernel_sec) && o.kernel_sec > 0 && o.raw_N_estimate > 0
            T(c).mass_ok = o.est_prior * N / (o.raw_N_estimate * max(1, o.kernel_sec * o.sampling_rate_used)) <= 5;
        end
    end
    % Validity, in order of preference: converged EM and both structural
    % checks (signal right of noise, signal mass consistent with the
    % detections), then just signal right of noise, then anything that fit.
    valid = arrayfun(@(t) isempty(t.error) && isfinite(t.bic) && t.signal_right_of_noise ...
                          && t.em_converged && t.mass_ok, T);
    if ~any(valid)
        valid = arrayfun(@(t) isempty(t.error) && isfinite(t.bic) && t.signal_right_of_noise, T);
    end
    if ~any(valid)
        valid = arrayfun(@(t) isempty(t.error) && isfinite(t.bic), T);
    end
    if ~any(valid), error('detection_stats:modelSelectFailed', 'No candidate produced a usable fit for component %d.', k); end
    tailScore = arrayfun(@(t) abs(log(max(t.tail_ratio_999, eps))), T);
    switch criterion
        case 'tail'
            score = tailScore;
        case 'ks'
            score = arrayfun(@(t) t.ks_noise, T);
            score(~isfinite(score)) = arrayfun(@(t) t.bic, T(~isfinite(score)));  % joint-EM branch has no KS
        otherwise   % 'bic' | 'aic': within a local-scale value; across values by tail calibration
            score = arrayfun(@(t) t.(criterion), T);
    end
    score(~valid) = Inf;
    spanKey = cellfun(@(v) sprintf('%g', double(v)), {C.local_scale_sec}, 'UniformOutput', false);
    ukeys = unique(spanKey);
    if numel(ukeys) > 1 && any(strcmp(criterion, {'bic','aic'}))
        % Likelihoods are not comparable across spans (different data), so
        % rank the spans by the family-free CV score of the block levels
        % (Inf/0 = global), then pick the family within the winning span.
        fs_k = hos(k).sampling_rate;
        if isempty(cache.scale_src), src_k = cache.xfilts(:,k); else, src_k = cache.scale_src(:,k); end
        [blk, ~, nBlk] = block_medians(src_k, cache.xthrs(:,k), fs_k, options.local_scale_block_sec);
        spanScore = inf(1, numel(ukeys)); best = nan(1, numel(ukeys));
        for u = 1:numel(ukeys)
            idx = find(strcmp(spanKey, ukeys{u}) & valid);
            if isempty(idx), continue, end
            sv = C(idx(1)).local_scale_sec;
            if ischar(sv), sv = res{idx(1)}.outs(k).local_scale.sec; end
            spanScore(u) = cv_span_score(blk, nBlk, options.local_scale_block_sec, sv);
            [~, j] = min(score(idx)); best(u) = idx(j);
            for i = idx, T(i).cv_span_score = spanScore(u); end
        end
        [~, u] = min(spanScore); sel = best(u);
    else
        [~, sel] = min(score);
    end
    T(sel).selected = true;
    o = res{sel}.outs(k);
    o.model_selection = T;
    % Effective parameters of the winner (model_select off, numeric span) so
    % that detection_stats(hos, x, outs(k).options) reproduces it directly;
    % the selection request itself is kept in options.requested.
    o.options.noise_dist = C(sel).noise_dist; o.options.signal_dist = C(sel).signal_dist;
    o.options.freeze_noise = C(sel).freeze_noise;
    o.options.local_scale_sec = o.local_scale.sec;
    o.options.model_select = false; o.options.select_criterion = criterion;
    o.options.requested = struct('local_scale_sec', options.local_scale_sec, 'model_select', model_select);
    if isempty(outs), outs = o; else, outs(k) = o; end
    pprobs(:,k) = res{sel}.pprobs(:,k);
    pprobs_outlier(:,k) = res{sel}.pprobs_outlier(:,k);
    local_scales(:,k) = res{sel}.local_scales(:,k);
    fprintf('  component %d: selected %s noise / %s signal, freeze=%d, local_scale_sec=%s (BIC %.1f, tail ratio %.2f)\n', ...
        k, C(sel).noise_dist, C(sel).signal_dist, C(sel).freeze_noise, num2str(o.local_scale.sec), ...
        T(sel).bic, T(sel).tail_ratio_999);
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
%   noise_dist  'gamma' | 'lognormal' | 'weibull' | 'chi2' | 'gengamma'
%   opts        optimset struct for fminsearch (gamma + Weibull)
%
% OUTPUTS
%   a_n, b_n            gamma shape, scale (NaN unless noise_dist='gamma')
%   mu_n, sig_n         lognormal mu, sigma (NaN unless 'lognormal')
%   a_w, b_w            Weibull scale, shape (NaN unless 'weibull');
%                       for 'gengamma' these carry scale a and shape d
%   nu_c                chi2 DF (NaN unless 'chi2'); for 'gengamma' the
%                       power p (the slots are only read back by
%                       make_noise_params_local)
%   mode_x              mode of the marginal in xrsm units
%   F_X_bulk, f_X_bulk  function handles for the marginal CDF and PDF
a_n = NaN; b_n = NaN; mu_n = NaN; sig_n = NaN;
a_w = NaN; b_w = NaN; nu_c = NaN;
sw = max(sum(w_obs), eps);
% Gamma and Weibull are fitted by exact WEIGHTED maximum likelihood on
% the full sample (one-dimensional root finds on the weighted sufficient
% statistics; see weighted_gamma_mle / weighted_weibull_mle). This
% replaces the earlier gamfit/wblfit on a random weighted subsample,
% which made every fit - and therefore the posterior - depend on the
% RNG state (runs were not reproducible to ~5 % in adjusted N).
switch noise_dist
    case 'gamma'
        try
            [a_n, b_n] = weighted_gamma_mle(x_obs.^2, w_obs);
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
        % Scaled chi-square (Satterthwaite): Z = x^2 ~ s * chi2(nu) with
        % nu = 2 E[Z]^2 / Var[Z] and s = E[Z] / nu, i.e. an average of
        % ~nu/2 half-normal squares. (The former unscaled chi2(E[Z]) had
        % variance 2 where the normalised statistic has ~2/n_eff.)
        z = x_obs.^2;
        mu_z  = sum(w_obs .* z) / sw;
        var_z = max(sum(w_obs .* (z - mu_z).^2) / sw, eps);
        nu_c  = [max(2 * mu_z^2 / var_z, 1), mu_z / max(2 * mu_z^2 / var_z, 1)];   % [nu, scale s]
        % f_X(x) = 2x * chi2pdf(x^2/s, nu)/s ; x_mode = sqrt(s (nu - 2) + s) = sqrt(s (nu-1))
        mode_x = sqrt(max(nu_c(2) * (nu_c(1) - 1), 0));
        F_X_bulk = @(x) chi2cdf(x.^2 / nu_c(2), nu_c(1));
        f_X_bulk = @(x) 2 .* x .* chi2pdf(x.^2 / nu_c(2), nu_c(1)) / nu_c(2);
    case 'weibull'
        try
            [a_w, b_w] = weighted_weibull_mle(x_obs, w_obs);
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
    case 'gengamma'
        % Stacy generalized gamma on x: scale a_g, shape d_g, power p_g.
        % Reported through the Weibull slots (a_w = a, b_w = d) and the
        % chi2 slot (nu_c = p); see make_noise_params_local.
        try
            [a_g, d_g, p_g] = weighted_gengamma_mle(x_obs, w_obs);
        catch
            a_g = max(sum(w_obs .* x_obs) / sw, 1e-3); d_g = 2; p_g = 2;   % Rayleigh-like fallback
        end
        a_w = a_g; b_w = d_g; nu_c = p_g;
        if d_g > 1
            mode_x = a_g * ((d_g - 1) / p_g)^(1 / p_g);
        else
            mode_x = 0;
        end
        F_X_bulk = @(x) gammainc((max(x, 0) / a_g).^p_g, d_g / p_g);
        f_X_bulk = @(x) exp(log(p_g) - d_g * log(a_g) - gammaln(d_g / p_g) ...
                            + (d_g - 1) .* log(max(x, realmin)) - (max(x, 0) / a_g).^p_g);
end
end

% =========================================================================
function [prm, it, converged, logL] = frozen_signal_em(family, xe, we, log_f_n, prm, maxit)
%FROZEN_SIGNAL_EM  EM over one parametric signal component (gamma on x^2 or
% lognormal on x) with the noise log-density log_f_n frozen, on (binned)
% samples xe with weights we. Closed-form weighted MLE M-steps.
W_total = sum(we); logL_prev = -Inf; logL = -Inf; it = 0;
lx = log(max(xe, eps));
for it = 1:maxit
    switch family
        case 'gamma'
            log_f_s = log(max(2 * xe .* gampdf(xe.^2, prm.a, prm.b), eps));
        otherwise
            log_f_s = -lx - log(prm.sigma) - 0.5 * log(2 * pi) - (lx - prm.mu).^2 / (2 * prm.sigma^2);
    end
    log_pi_n = log(max(1 - prm.pi_s, eps));
    log_pi_s = log(max(prm.pi_s, eps));
    m = max(log_pi_n + log_f_n, log_pi_s + log_f_s);
    lden = m + log(exp(log_pi_n + log_f_n - m) + exp(log_pi_s + log_f_s - m));
    gamma_t = exp(log_pi_s + log_f_s - lden);
    logL = sum(we .* lden);
    W_s = sum(we .* gamma_t);
    prm.pi_s = W_s / W_total;
    if W_s > 1
        switch family
            case 'gamma'
                z_eff   = xe.^2;
                mu_z    = sum(we .* gamma_t .* z_eff) / W_s;
                mlog_z  = sum(we .* gamma_t .* log(max(z_eff, eps))) / W_s;
                s_ = log(mu_z) - mlog_z;
                if s_ > 1e-10
                    a_s = (3 - s_ + sqrt((s_ - 3)^2 + 24 * s_)) / (12 * s_);
                    for ni = 1:50
                        f  = log(a_s) - psi(a_s) - s_;
                        fp = 1 / a_s - psi(1, a_s);
                        a_new = a_s - f / fp;
                        if abs(a_new - a_s) < 1e-6 * abs(a_s) + 1e-9, a_s = a_new; break, end
                        a_s = max(a_new, 1e-3);
                    end
                    prm.a = a_s; prm.b = mu_z / a_s;
                end
            otherwise
                prm.mu    = sum(we .* gamma_t .* lx) / W_s;
                prm.sigma = sqrt(max(sum(we .* gamma_t .* (lx - prm.mu).^2) / W_s, 1e-6));
        end
    end
    if it > 1 && abs(logL - logL_prev) < 1e-6 * max(abs(logL), 1), break, end
    logL_prev = logL;
end
converged = it < maxit;
end

function m = signal_mode_of(family, prm)
%SIGNAL_MODE_OF  Mode in x of the signal component.
switch family
    case 'gamma'
        m = sqrt(max(prm.b * (2 * prm.a - 1) / 2, 0));
    otherwise
        m = exp(prm.mu - prm.sigma^2);
end
end

function m = signal_median_of(family, prm)
%SIGNAL_MEDIAN_OF  Median in x of the signal component (location measure
% that, unlike the mode, is not pulled to zero by a wide lognormal).
switch family
    case 'gamma'
        m = sqrt(max(gaminv(0.5, prm.a, prm.b), 0));
    otherwise
        m = exp(prm.mu);
end
end

function q = noise_quantile(F, p, xmax)
%NOISE_QUANTILE  x at which the frozen noise CDF F reaches p (grid + interpolation).
xs = linspace(0, max(xmax, 1), 20000);
Fx = F(xs); Fx(~isfinite(Fx)) = 0;
j = find(Fx >= p, 1, 'first');
if isempty(j), q = xs(end); elseif j == 1, q = xs(1);
else, q = xs(j-1) + (p - Fx(j-1)) / max(Fx(j) - Fx(j-1), eps) * (xs(j) - xs(j-1));
end
end

% =========================================================================
function [a, b] = weighted_gamma_mle(z, w)
%WEIGHTED_GAMMA_MLE  Gamma(shape a, scale b) MLE with sample weights w >= 0.
% Score equations: log(a) - psi(a) = log(mean_w z) - mean_w(log z) =: s,
% b = mean_w(z) / a. The left side is monotone in a, so a single fzero on
% log(a) suffices; Minka's closed-form start is used as the bracket centre.
keep = isfinite(z) & z > 0 & isfinite(w) & w > 0;
z = z(keep); w = w(keep);
if numel(z) < 10, error('weighted_gamma_mle:tooFew', 'too few samples'); end
sw = sum(w);
m  = sum(w .* z) / sw;
ml = sum(w .* log(z)) / sw;
s  = log(m) - ml;
if ~(s > 1e-10)            % essentially degenerate (all equal): Gaussian limit
    a = 1e6; b = m / a; return
end
a0 = (3 - s + sqrt((s - 3)^2 + 24 * s)) / (12 * s);
f  = @(u) u - psi(exp(u)) - s;      % u = log(a)
lo = log(a0) - 2; hi = log(a0) + 2;
while f(lo) < 0, lo = lo - 2; if lo < -20, break, end, end
while f(hi) > 0, hi = hi + 2; if hi > 25, break, end, end
try
    u = fzero(f, [lo, hi]);
catch
    u = log(a0);
end
a = exp(u);
b = m / a;
end

function [lambda, k] = weighted_weibull_mle(x, w)
%WEIGHTED_WEIBULL_MLE  Weibull(scale lambda, shape k) MLE with weights.
% Profile score for k:
%   sum(w x^k log x)/sum(w x^k) - 1/k - mean_w(log x) = 0,
% lambda = (mean_w(x^k))^(1/k). Data are rescaled by their weighted
% geometric mean so that x^k cannot overflow.
keep = isfinite(x) & x > 0 & isfinite(w) & w > 0;
x = x(keep); w = w(keep);
if numel(x) < 10, error('weighted_weibull_mle:tooFew', 'too few samples'); end
sw = sum(w);
lx = log(x);
ml = sum(w .* lx) / sw;
c  = exp(ml);                       % weighted geometric mean
y  = x / c; ly = lx - ml;           % mean_w(ly) = 0
sd = sqrt(max(sum(w .* ly.^2) / sw, 1e-12));
k0 = 1.2825 / sd;                   % SD(log X) = pi/(k sqrt 6) for a Weibull
g  = @(k) sum(w .* y.^k .* ly) / max(sum(w .* y.^k), realmin) - 1 / k;   % increasing in k
lo = k0 / 4; hi = k0 * 4;
while g(lo) > 0 && lo > 1e-3, lo = lo / 2; end
while g(hi) < 0 && hi < 1e3,  hi = hi * 2; end
try
    k = fzero(g, [lo, hi]);
catch
    k = k0;
end
lambda = c * (sum(w .* y.^k) / sw)^(1 / k);
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
    case 'gengamma'
        s.a = a_wbl; s.d = b_wbl; s.p = nu_chi2(1);
        % d -> Inf with p -> 0 is the lognormal limit of the family: the
        % likelihood is flat along that ridge, so a fit that ran to the
        % power bound is "a lognormal by another name" and its a/d/p are
        % not individually meaningful. Say so.
        if s.p <= 0.05 && s.d >= 100
            s.limit = 'lognormal (p at lower bound, d large): equivalent to noise_dist=lognormal';
        elseif abs(s.p - 1) < 0.05
            s.limit = 'near gamma-on-x (p ~ 1)';
        elseif abs(s.d - s.p) < 0.05 * max(s.p, 1)
            s.limit = 'near Weibull (d ~ p)';
        end
    case 'chi2'
        s.nu = nu_chi2(1);
        if numel(nu_chi2) > 1, s.scale = nu_chi2(2); else, s.scale = 1; end
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
