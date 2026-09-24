function varargout = zerolag_bandcorr(hos,varargin)

% [R,fc,out] = zerolag_bandcorr(hos,[value],...)
%
% Zero-lag correlation of the rows of the Wigner-Ville distribution at
% carriers fc with each other (a 4th-order hosobject) or with the signal
% itself (a 3rd-order hosobject), summed from the reduced coefficients with
% the matrices ZEROLAG_INTEGRATOR builds:
%
%    order 4:  R(f1,f2) = sum_nu T(f1+nu/2, -f1+nu/2, -f2-nu/2, f2-nu/2)
%    order 3:  R(f)     = sum_nu B(f+nu/2, -f+nu/2)
%
% Inputs:
%   hos    a hosobject fitted with diagonal_slice = false. An array (a
%          chain of components) gives R(:,:,k) (order 4) or R(:,k) (order
%          3) per component; the matrices are rebuilt only when a
%          component's coefficient vector has a different length.
%   value  which coefficients to sum: 'bicoh' (default; the bias-corrected
%          polycoherence hos.bicohreduced, what MODGRAM plots), 'B' (the
%          raw moment spectrum hos.B), or a numeric vector the length of
%          hos.B (e.g. a significance-masked copy; a matrix with one column
%          per component is accepted).
%   'Normalize'  'sum' (default; the plain sum over the terms present),
%          'mean' (R./count, the mean coefficient along the line, NaN where
%          nothing was summed) or 'corr' (order 4, value 'B': R./sqrt(d*d')
%          with d the diagonal of the SUMS, the Gram matrix of the rows, so
%          the result is a correlation coefficient in [-1,1]; the diagonal
%          cells hold the pairs w1+w3 = 0 that the object's shighpass
%          removes, so this needs a fit with shighpass below
%          -0.25*sampling_rate/buffersize, and 'PairMean' to be a
%          covariance rather than a moment, see below).
%   'PairMean'  the struct ZEROLAG_PAIRMEAN returns (or a cell, one per
%          component): subtracts the second-order term
%          sum_nu E[M_f1(nu)]E[M_f2(nu)]* (order 4) or E[M_f(nu)]E[X(nu)]*
%          (order 3), which under a data window makes independent bands
%          look co-modulated. Meaningful for value 'B'.
%   Any other name-value pair ('Grid', 'MaxModulation', 'MinModulation')
%   is passed to ZEROLAG_INTEGRATOR.
%
% Outputs:
%   R      n x n x ncomp (order 4) or n x ncomp (order 3); real up to the
%          object's own symmetry (see ZEROLAG_INTEGRATOR on asymmetric
%          masks). Cells nothing was summed into are NaN.
%   fc     carrier axis (Hz).
%   out    struct: count (the terms summed per cell), I, Iconj (the
%          matrices), info (ZEROLAG_INTEGRATOR's), value, normalize.
%
% With no output argument, plots real(R) per component (an image at order
% 4, a line at order 3, empty cells transparent), in the manner of MODGRAM.
%
% Example, a batch estimate of the band-covariance matrix from a signal x:
%
%    hos = hosobject(4,N,fs,lowpass,[],[],'slowpass',numax,'diagonal_slice',false,'shighpass',-fs/N);
%    idx = (1:N)' + (0:floor((numel(x)-N)/(N/2)))*N/2;   % 50% overlap
%    FX = fft(x(idx).*hos.win);
%    hos.update_bispectrum(FX,true);
%    [~,~,info] = zerolag_integrator(hos);
%    S = zerolag_pairmean(FX,hos,info);
%    C = zerolag_bandcorr(hos,'B','Normalize','corr','PairMean',S,'Grid','full');
%
% See also ZEROLAG_INTEGRATOR, ZEROLAG_PAIRMEAN, HOSOBJECT/MODGRAM.

% C. Kovach 2026

value = 'bicoh';
if ~isempty(varargin) && (isnumeric(varargin{1}) || ...
        ((ischar(varargin{1}) || isstring(varargin{1})) && ismember(lower(char(varargin{1})),{'b','bicoh'})))
    value = varargin{1};
    varargin(1) = [];
end
p = inputParser;
p.FunctionName = mfilename;
p.KeepUnmatched = true;
p.addParameter('Normalize','sum',@(x)ischar(x)||isstring(x));
p.addParameter('PairMean',[],@(x)isempty(x)||isstruct(x)||iscell(x));
p.parse(varargin{:});
normalize = lower(char(p.Results.Normalize));
pairmean = p.Results.PairMean;
if ~ismember(normalize,{'sum','mean','corr'})
    error('%s: ''Normalize'' must be ''sum'', ''mean'' or ''corr'', not ''%s''.',mfilename,normalize);
end
intargs = [fieldnames(p.Unmatched)'; struct2cell(p.Unmatched)'];
intargs = intargs(:)';

order = double(hos(1).order);
if strcmp(normalize,'corr') && order == 3
    error('%s: ''corr'' needs the rows'' norms, which a 3rd-order object does not hold; use ''sum'' or ''mean''.',mfilename);
end

ncomp = numel(hos);
R = [];
lenprev = -1;
I = []; Iconj = []; info = [];
for k = 1:ncomp
    if isnumeric(value)
        v = value(:,min(k,size(value,2)));
        v = v(:);
    else
        switch lower(char(value))
            case 'b'
                v = hos(k).B(:);
            case 'bicoh'
                v = hos(k).bicohreduced(:);
        end
    end
    lenk = numel(v)-1;
    if lenk ~= lenprev
        [I,Iconj,info] = zerolag_integrator(hos(k),intargs{:});
        if info.len ~= lenk
            error('%s: the value vector has %d entries but hos(%d) holds %d reduced coefficients (+1).',...
                mfilename,numel(v),k,info.len);
        end
        lenprev = lenk;
    end
    v = v(1:end-1);
    v(isnan(v)) = 0;
    Rk = I*v + Iconj*conj(v);
    if ~isempty(pairmean)
        if iscell(pairmean)
            S = pairmean{min(k,numel(pairmean))};
        else
            S = pairmean;
        end
        Rk = Rk - secondorder(S,info);
    end
    if order == 4
        Rk = reshape(Rk,info.n,info.n);
    end
    switch normalize
        case 'mean'
            Rk = Rk./info.count;
        case 'corr'
            d = real(diag(Rk));
            Rk = Rk./sqrt(d*d');
    end
    Rk(info.count==0) = nan;
    if order == 4
        R = cat(3,R,Rk);
    else
        R = cat(2,R,Rk(:));
    end
end
fc = info.fc;
out.count = info.count;
out.I = I;
out.Iconj = Iconj;
out.info = info;
out.value = value;
out.normalize = normalize;

if nargout == 0
    if isnumeric(value)
        vlabel = 'input';
    else
        vlabel = char(value);
    end
    for k = 1:ncomp
        ax = subplot(1,ncomp,k);
        if order == 4
            Rp = real(R(:,:,k));
            imh = imagesc(ax,fc,fc,Rp);
            set(imh,'AlphaData',~isnan(Rp));
            axis(ax,'xy','square');
            colorbar(ax);
            xlabel(ax,'Band freq. (Hz)');
            ylabel(ax,'Band freq. (Hz)');
            ttl = 'Zero-lag WVD band correlation';
        else
            plot(ax,fc,real(R(:,k)));
            xlabel(ax,'Band freq. (Hz)');
            ylabel(ax,sprintf('%s of %s',normalize,vlabel));
            ttl = 'Zero-lag WVD-to-signal correlation';
        end
        if k == 1
            title(ax,sprintf('%s (%s of %s)',ttl,normalize,vlabel));
        else
            title(ax,sprintf('Residual after comp. %d',k-1));
        end
    end
else
    varargout = {R,fc,out};
end
end

function R2 = secondorder(S,info)
% sum_nu present(cell,nu) * E[M_f1(nu)] E[M_f2(nu)]*  (order 4)
% sum_nu present(f,nu)    * E[M_f(nu)]  E[X(nu)]*     (order 3), on the half grid, then folded
if ~isequal(size(S.m),[numel(info.fchalf) numel(info.nu)])
    error('zerolag_bandcorr:pairMean','''PairMean'' is %s but the lattice is %d carriers x %d modulation frequencies; build it with zerolag_pairmean(FX,hos,info) from the same info.',...
        mat2str(size(S.m)),numel(info.fchalf),numel(info.nu));
end
m = S.m;
m(isnan(m)) = 0;
nh = numel(info.fchalf);
pres = info.presenthalf;
switch info.order
    case 4
        R2 = zeros(nh^2,1);
        for q = 1:numel(info.nu)
            R2 = R2 + pres(:,q).*reshape(m(:,q)*m(:,q)',[],1);
        end
    case 3
        R2 = sum(pres.*(m.*conj(S.xbar(:)')),2);
end
if ~isempty(info.fold)
    R2 = info.fold*R2;
end
R2 = full(R2);
end
