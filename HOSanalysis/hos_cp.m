function out = hos_cp(hos, R, varargin)

% out = hos_cp(hos, R, ...)
%
% Rank-R CANDECOMP / PARAFAC (CP) decomposition of the higher-order
% spectra accumulated in the hosobject, hos (or of a signal given with
% 'X'), regarded as a decomposition of a full K-way
% moment tensor where all coefficients outside the subdomain of the  time-invariant 
% slice specified in the hosobject, hos, are missing entries.
% Complex factors are fitted by alternating least squares (ALS) over the
% observed entries only. 
%
%
% Input:
%
%   hos  - hosobject (the first element of an array is used) giving the
%          order, the sampling rate, the frequency indexing and, once it
%          has been fitted, the moment itself: its B, the average over its
%          tapered windows of its input (windows holding NaN samples
%          dropped). For a fitted hosd, hos(1) holds the moment of the
%          signal the fit started from.
%   R    - rank (number of components).
%
% Options (name, value):
%
%   'X'             - a signal: estimate the moment from it instead, in one
%                     batch through a reset copy of hos (get_input: its
%                     windows, NaN windows dropped). Needed when hos has not
%                     been fitted; for a fitted hos, its own input gives its
%                     B back exactly.
%   'Taper'         - with 'X': apply the hosobject's window, as its fit
%                     does (default true). Rectangular windows cut the
%                     signal at their edges, and the discontinuities spread
%                     over the whole subdomain.
%   'B'             - the reduced moment, numel(hos.freqindx.Is(:,1))
%                     entries (or one more: hosobject's trailing entry is
%                     ignored), instead of the hosobject's.
%   'Weights'       - one finite, non-negative weight per reduced
%                     coefficient (indexed like B: numel(hos.freqindx.Is(:,1))
%                     entries, or one more, the trailing one ignored) for a
%                     WEIGHTED least-squares fit. Every entry expanded from a
%                     reduced coefficient carries its weight, in the row-wise
%                     normal equations, G = sum w conj(W) W and rhs =
%                     sum w conj(W) V, and in the residual (see 'residual').
%                     The weights are scaled to mean 1 inside (the solution
%                     does not depend on their scale). [] (default): every
%                     entry weight 1, the unweighted fit, unchanged. Two
%                     natural choices: 1 ./ D.^2, with D = hos.D(1:n) the
%                     normalization of hos's polycoherence, fits the
%                     polycoherence B ./ D (the quantity HOSD works on); the
%                     inverse variance over the fit's windows of the
%                     per-window products that average to B fits every entry
%                     at its own noise level. Not to be confused with the
%                     output's 'weights', the component weights.
%   'Init'          - [] (default): random complex starts; an N x R0 matrix
%                     of waveforms on the hosobject's sample grid (the order
%                     of hos.feature: sample 0 first); or a hosobject array,
%                     whose .feature waveforms are taken. Each waveform's
%                     spectrum starts the corresponding column of every
%                     factor; columns beyond R0 start at random.
%   'NumStarts'     - random starts (default 8; 1 when 'Init' is given, in
%                     which case further starts are random). The start with
%                     the lowest residual is kept.
%   'MaxIterations' - ALS sweeps per start (default 4000).
%   'Tol'           - stop when the relative residual changes by less than
%                     Tol (relative) between checks (default 1e-6).
%   'CheckEvery'    - sweeps between residual checks (default 10).
%   'Seed'          - random seed (default 1; [] uses the global stream).
%                     The same seed gives the same result.
%   'Verbose'       - print one line per start (default false).
%
% Output: struct with fields
%
%   order, rank      K and R
%   freqs            the signed frequency bins of the factors (Hz), nb x 1
%   factors          1 x K cell of nb x R complex factors (the first K-1
%                    columns unit norm, the norms pushed into the last)
%   weights          1 x R, the norm of each column of the last factor
%   waveforms        nfft x R real waveforms on the hosobject's sample grid
%                    (sample 0 first, as hos.feature), from the first factor:
%                    inverse DFT, rotated by the global phase that makes it
%                    most nearly real (1/2 arg sum z^2), scaled by
%                    weight^(1/K). The fit does not identify a time shift
%                    (each factor times exp(i 2 pi f tau) leaves every entry
%                    unchanged): a seeded column is shifted (and its sign
%                    chosen) to match its seed's waveform, a random one is
%                    centered on its energy centroid.
%   residual         relative residual over the observed entries,
%                    sqrt(sum w |V - M|^2 / sum w |V|^2) (every w 1 without
%                    'Weights')
%   residualUnweighted  the same with every w 1 (equal to residual without
%                    'Weights'): the plain misfit of a weighted fit
%   weighted         true when 'Weights' was given
%   entryWeights     the weights as given, one per reduced coefficient ([]
%                    without 'Weights')
%   iterations       ALS sweeps run
%   settle           first checked sweep at which the residual is within
%                    0.1 % of its final value (a speed measure independent of
%                    the stopping rule)
%   history          rows [sweep, residual] at each check
%   starts           per start: init ('seeds' | 'random'), residual,
%                    iterations, settle, time, waveforms
%   startUsed        the start kept
%   nEntries         observed entries after expansion; nReduced, the
%                    reduced coefficients they came from
%   B                the reduced moment fitted
%   source           where it came from: "fit" (the hosobject's own B),
%                    "signal" ('X') or "given" ('B')
%
% Notes: the factors are not constrained to the conjugate symmetry of a real
% waveform's spectrum; multiplying each factor by exp(z f) also leaves every
% entry unchanged, so a real z (a spectral tilt) is not identified either
% and random starts land at small, start-dependent tilts. Started from HOSD's
% waveforms the fit settles much sooner than from random starts and refines
% them; it can also move a component away from a better HOSD estimate.
% The ALS updates every row of a factor at once: an R x R system per
% frequency bin (pagemldivide), summed over that bin's own observed entries
% through a sparse selection matrix.
%
% See also hosd/component_hos, hosobject/get_input, zerolag_integrator

% C. Kovach 2026

p = inputParser;
p.addParameter('X', []);
p.addParameter('B', []);
p.addParameter('Weights', []);
p.addParameter('Taper', true);
p.addParameter('Init', []);
p.addParameter('NumStarts', []);
p.addParameter('MaxIterations', 4000);
p.addParameter('Tol', 1e-6);
p.addParameter('CheckEvery', 10);
p.addParameter('Seed', 1);
p.addParameter('Verbose', false);
p.parse(varargin{:});
opt = p.Results;

h = hos(1);
K = double(h.order);
fx = h.freqindx;
Is = double(fx.Is);
nIs = size(Is, 1);
fs = double(h.sampling_rate);
fr = cellfun(@(f) double(f(:)), h.freqs, 'UniformOutput', false);
nfft = numel(fr{1});
df = fs / nfft;

%%% The reduced moment: given, re-estimated from a signal, or the one the
%%% hosobject accumulated while it was fitted (window_number counts its
%%% windows; an unfitted hosobject's B is not a moment)
if ~isempty(opt.B)
    B = double(opt.B(:));
    source = "given";
elseif ~isempty(opt.X)
    hh = hosobject(h);
    hh.reset;
    hh.get_input(double(opt.X(:)), logical(opt.Taper), false, true);
    B = double(gather(hh.B(:)));
    source = "signal";
elseif h.window_number > 0
    B = double(gather(h.B(:)));
    source = "fit";
else
    error('hos_cp:noMoment', ['The hosobject holds no moment (it has not been fitted): ' ...
        'give a signal (''X'') or the reduced moment (''B'').'])
end
if numel(B) < nIs
    error('hos_cp:badB', 'B has %d entries for %d reduced coefficients.', numel(B), nIs)
end
B = B(1:nIs);

%%% Entry weights: one per reduced coefficient, indexed like B ([] = unweighted)
wEntry = [];
if ~isempty(opt.Weights)
    wEntry = double(opt.Weights(:));
    if numel(wEntry) < nIs
        error('hos_cp:badWeights', 'Weights has %d entries for %d reduced coefficients.', numel(wEntry), nIs)
    end
    wEntry = wEntry(1:nIs);
    if any(~isfinite(wEntry)) || any(wEntry < 0)
        error('hos_cp:badWeights', 'Weights must be finite and non-negative.')
    end
end

%%% Observed entries: each reduced coefficient as a tensor index over signed
%%% bins, the last axis holding conj(X(sum)) = X(-sum); then every
%%% permutation and every conjugate entry
Mr = zeros(nIs, K);
for k = 1:K
    Mr(:, k) = round(fr{k}(Is(:, k)) / df);
end
Mr(:, K) = -Mr(:, K);
ok = isfinite(B);
Mr = Mr(ok, :);
vr = B(ok);
wr = [];
if ~isempty(wEntry)
    wr = wEntry(ok);
    if ~any(wr > 0)
        error('hos_cp:badWeights', 'Every weight of an observed entry is zero.')
    end
    wr = wr / mean(wr);
end
P = perms(1:K);
nP = size(P, 1);
n0 = size(Mr, 1);
E = zeros(2 * nP * n0, K);
V = complex(zeros(2 * nP * n0, 1));
Wt = zeros(2 * nP * n0, 1);
for q = 1:nP
    E((q - 1) * n0 + (1:n0), :) = Mr(:, P(q, :));
    V((q - 1) * n0 + (1:n0)) = vr;
    E((nP + q - 1) * n0 + (1:n0), :) = -Mr(:, P(q, :));
    V((nP + q - 1) * n0 + (1:n0)) = conj(vr);
    if ~isempty(wr)
        Wt((q - 1) * n0 + (1:n0)) = wr;
        Wt((nP + q - 1) * n0 + (1:n0)) = wr;
    end
end
[E, ia] = unique(E, 'rows', 'stable');
V = V(ia);
if isempty(wr)
    Wt = [];
else
    Wt = Wt(ia);
end
nobs = size(E, 1);
mMax = max(abs(E(:)));
nb = 2 * mMax + 1;
idx = E + mMax + 1;
S = cell(1, K);
for q = 1:K
    S{q} = sparse(idx(:, q), (1:nobs)', 1, nb, nobs);
end
nV = norm(V);
nVw = nV;
if ~isempty(Wt)
    nVw = sqrt(sum(Wt .* abs(V).^2));
end

%%% Starting points
if isa(opt.Init, 'hosobject')
    W0 = zeros(nfft, numel(opt.Init));
    for k = 1:numel(opt.Init)
        w = double(opt.Init(k).feature(:));
        W0(1:min(nfft, numel(w)), k) = w(1:min(nfft, numel(w)));
    end
else
    W0 = double(opt.Init);
end
if ~isempty(W0) && size(W0, 1) ~= nfft
    w = zeros(nfft, size(W0, 2));
    w(1:min(nfft, size(W0, 1)), :) = W0(1:min(nfft, size(W0, 1)), :);
    W0 = w;
end
W0 = W0(:, 1:min(R, size(W0, 2)));
nSeed = size(W0, 2);
if isempty(opt.NumStarts)
    if nSeed > 0, nStarts = 1; else, nStarts = 8; end
else
    nStarts = opt.NumStarts;
end
if isempty(opt.Seed)
    rs = RandStream.getGlobalStream;
else
    rs = RandStream('mt19937ar', 'Seed', opt.Seed);
end
mBins = (-mMax:mMax)';
A0 = [];
if nSeed > 0
    Aw = fft(W0);
    A0 = Aw(mod(mBins, nfft) + 1, :);
end

%%% Fit each start
starts = repmat(struct('init', "", 'residual', NaN, 'iterations', 0, 'settle', NaN, ...
    'time', NaN, 'waveforms', []), 1, nStarts);
best = Inf;
for st = 1:nStarts
    F = cell(1, K);
    for q = 1:K
        F{q} = complex(randn(rs, nb, R), randn(rs, nb, R));
    end
    if st == 1 && nSeed > 0
        for q = 1:K
            F{q}(:, 1:nSeed) = A0;
        end
        starts(st).init = "seeds";
    else
        starts(st).init = "random";
    end
    t0 = tic;
    [F, res, its, hist, resU] = als(F, idx, S, V, nV, nVw, Wt, R, K, nb, opt);
    starts(st).time = toc(t0);
    starts(st).residual = res;
    starts(st).iterations = its;
    starts(st).settle = settleIteration(hist);
    seeds = [];
    if starts(st).init == "seeds"; seeds = W0; end
    starts(st).waveforms = toWaveforms(F, mBins, nfft, K, seeds);
    if opt.Verbose
        fprintf('hos_cp start %d (%s): residual %.3g, %d sweeps, settles by %d, %.2f s\n', st, ...
            starts(st).init, res, its, starts(st).settle, starts(st).time);
    end
    if res < best
        best = res;
        out = struct('order', K, 'rank', R, 'freqs', mBins * df, 'factors', {F}, ...
            'weights', vecnorm(F{K}), 'waveforms', starts(st).waveforms, 'residual', res, ...
            'residualUnweighted', resU, 'weighted', ~isempty(Wt), 'entryWeights', wEntry, ...
            'iterations', its, 'settle', starts(st).settle, 'history', hist, 'starts', [], ...
            'startUsed', st, 'nEntries', nobs, 'nReduced', n0, 'B', B, 'source', source);
    end
end
out.starts = starts;

end

% -------------------------------------------------------------------------
function [F, res, it, hist, resU] = als(F, idx, S, V, nV, nVw, Wt, R, K, nb, opt)
% ALS over the observed entries: each factor row is the (weighted, when Wt
% is not empty) least-squares solution over that row's own entries, an
% R x R system per row. res is the residual the fit minimizes (weighted
% with Wt), resU the unweighted one at the same point.
prev = NaN;
res = NaN;
resU = NaN;
hist = zeros(0, 2);
it = 0;
for it = 1:opt.MaxIterations
    for q = 1:K
        oth = [1:q-1, q+1:K];
        W = F{oth(1)}(idx(:, oth(1)), :);
        for o = oth(2:end)
            W = W .* F{o}(idx(:, o), :);
        end
        cW = conj(W);
        if ~isempty(Wt)
            cW = cW .* Wt;
        end
        G = complex(zeros(R, R, nb));
        rhs = complex(zeros(R, 1, nb));
        for r = 1:R
            rhs(r, 1, :) = reshape(S{q} * (cW(:, r) .* V), 1, 1, nb);
            for s = 1:R
                G(r, s, :) = reshape(S{q} * (cW(:, r) .* W(:, s)), 1, 1, nb);
            end
        end
        G = G + 1e-12 * max(abs(G(:))) * eye(R);
        F{q} = reshape(pagemldivide(G, rhs), R, nb).';
    end
    for q = 1:K-1
        n = vecnorm(F{q});
        n(n == 0) = 1;
        F{q} = F{q} ./ n;
        F{K} = F{K} .* n;
    end
    if mod(it, opt.CheckEvery) == 0 || it == opt.MaxIterations
        Mo = F{1}(idx(:, 1), :);
        for q = 2:K
            Mo = Mo .* F{q}(idx(:, q), :);
        end
        d = V - sum(Mo, 2);
        resU = norm(d) / nV;
        if isempty(Wt)
            res = resU;
        else
            res = sqrt(sum(Wt .* abs(d).^2)) / nVw;
        end
        hist(end + 1, :) = [it, res]; %#ok<AGROW>
        if abs(prev - res) < opt.Tol * max(res, realmin)
            break
        end
        prev = res;
    end
end
end

% -------------------------------------------------------------------------
function n = settleIteration(hist)
% First checked sweep within 0.1 % of the final residual.
if isempty(hist)
    n = NaN;
    return
end
n = hist(find(hist(:, 2) <= hist(end, 2) * (1 + 1e-3), 1), 1);
end

% -------------------------------------------------------------------------
function Y = toWaveforms(F, mBins, nfft, K, seeds)
% Real waveforms from the first factor (see the help: phase rotation,
% scaling, shift and sign).
R = size(F{1}, 2);
Y = zeros(nfft, R);
w = vecnorm(F{K});
n = (0:nfft-1)';
for r = 1:R
    Z = zeros(nfft, 1);
    Z(mod(mBins, nfft) + 1) = F{1}(:, r) * w(r)^(1/K);
    z = ifft(Z) * nfft;
    s2 = sum(z.^2);
    th = 0;
    if abs(s2) > 1e-9 * sum(abs(z).^2)
        th = angle(s2) / 2;
    end
    y = real(z * exp(-1i * th));
    if ~isempty(seeds) && r <= size(seeds, 2) && any(seeds(:, r))
        c = real(ifft(fft(seeds(:, r)) .* conj(fft(y))));
        [~, k] = max(abs(c));
        y = circshift(y, k - 1);
        if c(k) < 0; y = -y; end
    else
        ph = angle(sum((y.^2) .* exp(1i * 2 * pi * n / nfft)));
        y = circshift(y, -round(ph / (2 * pi) * nfft));
    end
    Y(:, r) = y;
end
end
