function [P, N, C] = laguerreBasis(tt, ord, tau, fs, sides, continuity, circular)
%LAGUERREBASIS  Discrete Laguerre basis over event lags: causal, anticausal or both.
%
%   [P, N, C] = laguerreBasis(tt, ord, tau, fs, sides, continuity, circular)
%
%   tt          lag axis in seconds (column); the sample with tt == 0 is the
%               event sample. For the full-record design matrix pass
%               tt = (0:n-1)'/fs with circular = true: anticausal lags then
%               wrap to the end of the record (circular convolution).
%   ord, tau    Laguerre order (ord+1 functions per side) and time constant (s)
%   fs          sampling rate (Hz)
%   sides       'causal'      lags >= 0 (the event sample and after)
%               'anticausal'  lags <= -1 sample (strictly before the event)
%               'both'        causal and anticausal sets side by side
%   continuity  for sides = 'both' only:
%               'none'  the two halves are free (the kernel may jump at 0)
%               'C0'    both halves extrapolate to the same value at the
%                       event: sum(b) == sum(c) (every function has the same
%                       value phi[0] there)
%               'C1'    ... and the same slope there, so the first
%                       differences are equal straight across the event:
%                       k[1]-k[0] == k[0]-k[-1]
%               Imposed by reparametrising the 2(ord+1) coefficients [b; c]
%               (b causal, c anticausal) with the null space of the
%               constraint rows, N = null(C): P = [Pc Pa] * N, the fitted
%               coefficients are gamma and [b; c] = N * gamma.
%   circular    true for the record-length design, false (default) for a window
%
%   P   numel(tt) x q basis, q = ord+1 (one side), 2(ord+1) (both, free),
%       2(ord+1)-1 (C0) or 2(ord+1)-2 (C1)
%   N   (ord+1 or 2(ord+1)) x q map from the fitted coefficients back to the
%       unconstrained ones (identity when unconstrained)
%   C   constraint rows on [b; c] (empty when unconstrained)
%
%   The causal functions are the impulse responses of laguerreFilt with its
%   one-sample delay removed, so every function starts at the event sample
%   with the same value phi[0]; the anticausal set is the same sequence
%   mirrored ABOUT THE EVENT, lag -m carrying phi[m] for m >= 1, so both
%   halves are samples of one function of |lag| that meet at the event
%   (the mirror is not shifted by a sample, which would make the halves
%   meet half a sample before the event and turn a 'C0' constraint into a
%   flat step there). Used by model.designMtx / model.get_event_window for
%   the 'laguerre' (and 'forward_laguerre' / 'backward_laguerre') time
%   bases.
%
%   See also LAGUERREFILT, MODEL.

% C. Kovach 2026

if nargin < 5 || isempty(sides);      sides = 'causal';   end
if nargin < 6 || isempty(continuity); continuity = 'none'; end
if nargin < 7 || isempty(circular);   circular = false;   end

tt = tt(:);
L  = numel(tt);
nb = ord + 1;
if L < 2
    error('laguerreBasis:tooShort', 'The lag axis needs at least two samples.');
end
if circular
    i0 = 1;
else
    [~, i0] = min(abs(tt));
    if abs(tt(i0)) > 0.5 / fs
        error('laguerreBasis:noEventSample', ...
            'The lag axis must contain the event sample (lag 0); it spans [%g %g] s.', tt(1), tt(end));
    end
end

h = laguerreFilt(L + 1, ord, tau, fs);   % impulse responses; row 1 is the filter's one-sample delay
h = h(2:end, :);                          % phi[0 .. L-1]

Pc = zeros(L, nb);                        % causal: phi[m] at lag +m
m  = 0:(L - i0);
Pc(i0 + m, :) = h(1 + m, :);

if circular
    Pa = circshift(flipud(Pc), [1 0]);    % phi[m] at index L+1-m, i.e. lag -m wrapped to the end
    Pa(1, :) = 0;                         % lag 0 belongs to the causal set
else
    Pa = zeros(L, nb);                    % anticausal: phi[m] at lag -m, m >= 1
    m  = 1:(i0 - 1);
    Pa(i0 - m, :) = h(1 + m, :);
end

switch lower(char(sides))
    case {'causal', 'forward', 'after'}
        P = Pc; N = eye(nb); C = zeros(0, nb);
    case {'anticausal', 'backward', 'before'}
        P = Pa; N = eye(nb); C = zeros(0, nb);
    case {'both', 'two-sided', 'twosided', 'around'}
        switch lower(char(continuity))
            case {'none', '', 'off'}
                C = zeros(0, 2 * nb);
            case {'c0', 'value'}
                C = [h(1, :), -h(1, :)];                    % both halves reach phi[0]'*b == phi[0]'*c at the event
            case {'c1', 'slope'}
                d = h(2, :) - h(1, :);
                C = [h(1, :), -h(1, :); d, d];              % ... with equal slopes: k[1]-k[0] == k[0]-k[-1]
            otherwise
                error('laguerreBasis:badContinuity', ...
                    'continuity must be ''none'', ''C0'' or ''C1'' (got ''%s'').', continuity);
        end
        if isempty(C)
            N = eye(2 * nb);
        else
            % Unit-norm rows: null() decides the rank relative to the
            % largest singular value, and the slope row grows like
            % ord/alpha while the value row stays O(1), so without this a
            % small alpha (short tau) would silently drop the value
            % constraint.
            C = C ./ sqrt(sum(C.^2, 2));
            N = null(C);
            if isempty(N)
                error('laguerreBasis:overConstrained', ...
                    'The %s continuity constraint leaves no free coefficient at order %d; raise the order.', ...
                    upper(continuity), ord);
            end
        end
        P = [Pc, Pa] * N;
    otherwise
        error('laguerreBasis:badSides', ...
            'sides must be ''causal'', ''anticausal'' or ''both'' (got ''%s'').', sides);
end
