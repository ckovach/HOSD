function [lin,inside] = square_index(me,W)

% [lin,inside] = square_index(me,W)
%
% Linear indices into the square form (freqindx.remap and the arrays that
% share its shape) of the points whose first order-1 frequencies are the
% rows of W (Hz, one column per axis). Each frequency is matched to the
% nearest value of its axis (freqindx.Bfreqs{k}), as a nearest-neighbor
% interpolation over the full array would match it; a frequency outside the
% axis range gives lin = 0 and inside = false. The coefficient at such a
% point is then v(remap(lin)), conjugated where PDconj(lin) is set, for a
% reduced-form vector v -- without forming the full array (see MODGRAM,
% PLOTPOLYCOH).
%
% See also FREQ2INDEX, MODGRAM, PLOTPOLYCOH.

% C. Kovach 2026

fx = me.freqindx;
sz = size(fx.remap);
K = numel(sz);
if size(W,2) ~= K
    error('square_index:size','W must have %d columns (one per square-form axis), not %d.',K,size(W,2));
end
n = size(W,1);
subs = cell(1,K);
inside = true(n,1);
for k = 1:K
    [f,srt] = sort(fx.Bfreqs{k}(:));
    idx = interp1(f,(1:numel(f))',W(:,k),'nearest');    % NaN outside the axis range
    ok = ~isnan(idx);
    s = zeros(n,1);
    s(ok) = srt(idx(ok));
    inside = inside & ok;
    subs{k} = s;
end
lin = zeros(n,1);
subs = cellfun(@(s)s(inside),subs,'uniformoutput',false);
lin(inside) = sub2ind(sz,subs{:});
