function [outs,pprobs,xdets,xsnrs,xrsms,xfsds,xfilts,xthrs,gs,pprobs_outlier,local_scales,xmisses] = detection_stats(hos, x, varargin)

%[outs,pprob,xdet,xsnr,xrsm] = detection_stats(hos, x, 'noise_dist','gamma', ...)
% Estimates a posterior probability on events detected by xdetect(hos,x)
% under a 2-component mixture model.
% DETECTION_STATS is HOSOBJECT/XDETECT followed by FIT_EVENT_MM, which fits
% the mixture to detections already made (for instance one component of a
% chain); see help fit_event_mm.
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
%                  0 (default) reproduces the previous behavior: one
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
%                  and summarized in outs(k).local_scale (.candidates,
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
%                  reconstruction). A component whose detection filter
%                  cannot be applied to the residual is measured on its
%                  xfilt instead, with a warning, and its
%                  outs(k).options.local_scale_source says 'xfilt'; the
%                  other components keep the residual (all of them fall
%                  back when the reconstruction itself fails).
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
%                  run, on each component separately. Joint-EM candidates
%                  (freeze_noise false) are opt-in: their free 'signal'
%                  component tends to absorb the heteroskedastic bulk of
%                  the noise and wins the likelihood with a prior of 10-20 %
%                  and an adjusted N several times that of the frozen fits,
%                  so they should only be compared against each other.
%                  Candidates that error, whose EM did not converge, whose
%                  signal component lies left of the noise, or whose signal
%                  mass is > 5x what the detections can account for are
%                  excluded (in that order, as long as something remains),
%                  each for the component on which it happens. The table of
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
%   'smoothing_window'  the window xdetect smooths the thresholded output
%                  with before picking peaks (its 4th input): 'auto'
%                  (default) = the window derived from each component's
%                  feature; 0 = no smoothing; a scalar = the standard
%                  deviation of a Hann window in seconds; a vector = the
%                  window (non-negative); a matrix = one window per
%                  component (e.g. the gs of an earlier call). The mixture
%                  is fit on the statistic smoothed with the same window.
%   'mask'         samples to leave out (xdetect's 5th input): 'auto'
%                  (default) = the NaN samples of x; a logical array with
%                  one row per sample marks further samples, which are set
%                  to NaN (the residual the running scale is measured on
%                  gets the same gaps). The detection filter output is
%                  undefined over the filter's support around each such
%                  sample (xmisses, xdetect's 8th output); xdetect leaves
%                  those samples out of the noise SD and the cumulant
%                  threshold, and here the smoothed statistic is NaN
%                  wherever its window reaches one, so those samples take
%                  no part in the noise or signal fit, the likelihood, the
%                  running scale or the cumulant readouts, and their
%                  posterior is 0. A detection within reach of one has NaN
%                  in outs(k).snrs. Recorded in outs(k).options.mask as
%                  given: pass a sparse logical, or NaN the data and leave
%                  'auto', to keep options small.
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
%                       nu = E[Z], which has variance 2 where the normalized
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
%       flexible signal family from cannibalizing noise mass.
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
%   outs   - struct array (one element per hos channel) of fit stats. The
%       components are fit independently of each other: outs(k) and the
%       k-th columns of pprob, pprobs_outlier and local_scales are those of
%       FIT_EVENT_MM on component k alone (see help fit_event_mm, ONE
%       COMPONENT OF A CHAIN). Under model_select, when the winners of
%       different components come from branches whose outs have different
%       fields, outs carries the union of their fields, [] where a
%       component's winner has no such field. Key distributional
%       parameter fields are grouped into substructs:
%         outs(k).noise_params   .type plus family-specific fitted params
%                                e.g. {gamma:.a,.b,.DF_eff}
%                                     {lognormal:.mu,.sigma}
%                                     {weibull:.a,.b}  {chi2:.nu}
%                                and .mode (x-space mode of fitted density)
%         outs(k).sig_params     analogous, with .type in {gamma,
%                                noncentral_gamma, lognormal, empirical}
%                                noncentral_gamma additionally has .lambda
%         outs(k).outlier_params .type='none' (default) or 'gpd' with
%                                .u, .k, .sigma, .pi, .pi_max. With
%                                outlier_dist 'gpd' in the joint EM, 'none'
%                                says that component k's GPD tail fit fell
%                                back to the 2-component model (fewer than
%                                50 exceedances, or gpfit failing); each
%                                component decides for itself.
%         outs(k).empirical_pdf  [length(snrx) x 1] empirical density of
%                                the data the EM was fit on, on the same
%                                px0 bin edges as sigpdf/noisepdf. For
%                                visual goodness-of-fit comparison.
%       Per-family parameters that don't apply to the chosen distribution
%       are simply not present in the substruct (no NaN clutter).
%       outs(k).options        every input parameter as EFFECTIVELY used
%                                (an 'auto' span is replaced by the chosen
%                                span, model_select is off, the selected
%                                families are filled in, local_scale_source
%                                is 'xfilt' when component k's residual
%                                could not be formed); pass it back as
%                                'options' to reproduce the fit without
%                                repeating the search. .requested keeps
%                                the original local_scale_sec/model_select.
%       outs(k).local_scale    .source ('residual' | 'xfilt', or 'none'
%                                when local_scale_sec is 0),
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
%   xmisses - [N x K] logical, true where a component's detection filter
%       output is undefined (see 'mask'); xfilts and xthrs are 0 there.
%
%See also HOSOBJECT/XDETECT, FIT_EVENT_MM

%C. Kovach 2026

% Every option is checked, and every warning about the options given,
% before the detection runs (the 'options' form of fit_event_mm).
[options, reads_residual] = fit_event_mm('options', varargin{:});
[xdets,xsnrs,xrsms,xfsds,xfilts,xthrs,gs,xmisses] = xdetect(hos,x,options.pow,options.smoothing_window,options.mask);
det = struct('xdet', xdets, 'xsnr', xsnrs, 'xrsm', xrsms, 'xfsd', xfsds, ...
             'xfilt', xfilts, 'xthr', xthrs, 'g', gs, 'xmiss', xmisses);
% Signal the running scale is measured on (see 'local_scale_source'),
% formed only when a span reads it: 'auto', a span > 0, or several spans
% (model selection ranks them on it). A model_select struct's
% local_scale_sec replaces local_scale_sec as the spans tried, and model
% selection at a zero span never reads it. The residual costs one xrec
% plus one filter pass per component.
if reads_residual
    % An explicit mask: xdetect has checked its shape and set those samples
    % to NaN on its own copy of x. The same here, so that the residual has
    % the same gaps.
    if ~ischar(options.mask), x(logical(options.mask)) = NaN; end
    det.xfilt_residual = residual_filter_output(hos, x, xfilts);
end
[outs, pprobs, pprobs_outlier, local_scales] = fit_event_mm(hos, det, options);
end

function xsrc = residual_filter_output(hos, x, xfilts)
%RESIDUAL_FILTER_OUTPUT  Each component's detection filter applied to the
% residual after removing EVERY component's reconstructed signal
% (hosobject/xrec: thresholded filter output * feature waveform, LMSE
% scaled). A component whose filter cannot be applied to the residual gets
% a NaN column, with a warning; all of them do when the reconstruction
% fails (e.g. hos is not a hosobject). FIT_EVENT_MM then measures that
% component's running scale on its xfilt.
xsrc = nan(size(xfilts), 'like', xfilts);
try
    xr = xrec(hos, x(:));                     % N x K (deflated, all components)
    xr(isnan(xr)) = 0;
    xres = x(:) - sum(xr, 2);
catch ME
    warning('detection_stats:residualScale', ...
        'Could not form the residual filter output (%s); the running scale is measured on xfilt instead.', ME.message);
    return
end
for k = 1:numel(hos)
    try
        xf = apply_filter(hos(k), xres, false, false);
        xsrc(:,k) = xf(:);
    catch ME
        warning('detection_stats:residualScale', ...
            ['Could not form the residual filter output of component %d (%s); its running scale ', ...
             'is measured on xfilt instead.'], k, ME.message);
    end
end
end
