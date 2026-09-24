function [I,Iconj,info] = zerolag_integrator(hos,varargin)

% [I,Iconj,info] = zerolag_integrator(hos,...)
%
% Summation matrices that integrate the reduced (non-redundant) 3rd- or
% 4th-order coefficients of a hosobject along the modulation-frequency
% direction, giving the zero-lag correlation of the rows of the signal's
% Wigner-Ville distribution (WVD) with each other (order 4) or with the
% signal itself (order 3). The WVD row at carrier f has Fourier transform
% M_f(nu) = X(f+nu/2)X*(f-nu/2) along time, so by Parseval
%
%   order 4:  R(f1,f2) = sum_nu E[M_f1(nu) M_f2(nu)*]
%                      = sum_nu T(f1+nu/2, -f1+nu/2, -f2-nu/2, f2-nu/2)
%   order 3:  R(f)     = sum_nu E[M_f(nu) X(nu)*]
%                      = sum_nu B(f+nu/2, -f+nu/2)
%
% with |nu| <= numax, where T(w1,w2,w3,w4) = E[X(w1)X(w2)X(w3)X(w4)] and
% B(w1,w2) = E[X(w1)X(w2)X*(w1+w2)] are the moment spectra accumulated in
% hos.B: the expected zero-lag inner product of two WVD rows, or of a WVD
% row with the slow (|nu| <= numax) part of the signal (the bispectrum
% being the cross-spectrum between the WVD and the signal). Read the other
% way, a square-form point (w1,w2,...) belongs to the carrier f = (w1-w2)/2
% (and f2 = (w4-w3)/2) at nu = w1+w2, a linear change of coordinates on the
% polyspectral domain, so every square-form point lies on exactly one of
% the lines summed here.
%
% Usage, in the manner of the Imats built by INTEGRATOR:
%
%    [I,Iconj,info] = zerolag_integrator(hos);
%    B = hos.B(1:end-1);                          % or hos.bicohreduced(1:end-1)
%    R = I*B + Iconj*conj(B);                     % order 3: n x 1
%    R = reshape(R, info.n, info.n);              % order 4: n x n
%    Rmean = R./info.count;                       % mean coefficient per cell
%
% The rows index carriers (order 3: row i is fc(i)) or carrier pairs (order
% 4: row = sub2ind([n n],i,j) for (f1,f2) = (fc(i),fc(j))); the columns
% index the reduced vector hos.B(1:end-1). Coefficients the object never
% estimated (masked out by lowpass / highpass / slowpass / shighpass /
% xlowpass / hos.mask, i.e. those hos.freqindx.remap sends to the trailing
% NaN slot) are simply absent, so info.count, the sum of the weights in
% each row, says how many terms each cell actually sums. A square-form
% point that is the conjugate of its principal-domain representative is
% summed through Iconj (applied to conj(B)), the others through I; one
% reduced coefficient stands for a symmetry class of square-form points and
% contributes to every cell whose line passes through one of them, with the
% multiplicity recorded in the matrices. Under the object's own
% (sign-symmetric) masks the nu and -nu terms of a cell are conjugates, so
% I and Iconj carry equal counts and R is real; an asymmetric hos.mask shows
% up as an imaginary part.
%
% On the DFT grid (spacing df = sampling_rate/fftN) a product X(k1)X*(k2)
% has carrier (k1+k2)df/2 and modulation (k1-k2)df: the discrete WVD lives
% on a carrier grid of spacing df/2, carriers ON the DFT grid see the even
% multiples of df as nu and carriers midway between bins the odd multiples,
% and (order 4) a cell (f1,f2) is populated only when f1 and f2 lie on the
% same lattice. Options:
%
%   'Grid'   'half' (default): carriers fc = 0:df/2:lowpass, every
%            coefficient with weight 1. Order-4 cells off the two lattices
%            (f1 on the DFT grid and f2 midway, or the reverse) are empty
%            (count 0). This is the exact summation matrix.
%            'full': carriers on the DFT grid, fc = 0:df:lowpass. Even-nu
%            coefficients keep weight 1 at their own cell; an odd-nu
%            coefficient, whose carrier sits midway between DFT bins, is
%            split over the neighboring cells (weight 1/2 to each of the
%            two neighbors at order 3, 1/4 to each of the four at order 4:
%            a linear / bilinear fold of the half grid). Every cell is
%            populated and the whole tube is used, at the price of a
%            +/- df/2 smoothing of the carrier axes.
%   'MaxModulation'  largest |nu| summed, in Hz. Default: the object's own
%            limit, so that every cell sums the same set of nu (up to the
%            edges and the highpass): slowpass - 0.25*sampling_rate/buffersize
%            at order 4 (the pair-sum lowpass with FIND_PRINCIPAL_DOMAIN's
%            tolerance; cells within slowpass of the diagonal would
%            otherwise hold extra terms through the (13)(24) pairing), and
%            just below xlowpass at order 3 (the "OR" lowpass, |w3| = |nu| <
%            xlowpass), or lowpass when no xlowpass is set.
%   'MinModulation'  smallest |nu| summed, in Hz (default 0; the object's
%            shighpass (order 4) or highpass on w3 (order 3) already removes
%            the small pair sums, nu = 0 included, which carry the
%            power-spectrum product).
%
% info fields: fc (carrier axis, Hz), n, count (n x n or n x 1), nu (the
% modulation frequencies enumerated, Hz), nupresent (those with at least
% one estimated coefficient), present (sparse rows x
% numel(nu): the weight of each (cell,nu) term) and presenthalf (the same
% on the half grid, for ZEROLAG_BANDCORR's second-order correction), fold
% (the 'full' fold matrix, [] for 'half'), fchalf (the half-grid
% carriers), grid, order, len (number of reduced coefficients), tube
% (number of distinct reduced coefficients used), df.
%
% The object must have been fitted with diagonal_slice = false: the
% diagonal slice meets the line of cell (f1,f2) only at the single point
% nu = f2-f1. With it off, the slowpass tube an order-4 object retains
% (every frequency has a partner with |w_i + w_j| <= slowpass) is exactly
% the set of coefficients these sums use, O(Nf^2 Nnu) of them against Nf^3
% for the full trispectrum and Nf Nnu for the modulogram (the diagonal
% slice of the same tube); at order 3 the xlowpass region (some |w_k| <
% xlowpass) plays the same role, O(Nf Nnu) coefficients. The matrices are
% built in one vectorized pass over the (f1,f2,nu) lattice and applied as a
% sparse product, both O(nnz).
%
% See also INTEGRATOR, ZEROLAG_BANDCORR, ZEROLAG_PAIRMEAN, HOSOBJECT/MODGRAM,
% FREQ2INDEX, FIND_PRINCIPAL_DOMAIN.

% C. Kovach 2026

p = inputParser;
p.FunctionName = mfilename;
p.addParameter('Grid','half',@(x)ischar(x)||isstring(x));
p.addParameter('MaxModulation',[],@(x)isempty(x)||(isscalar(x)&&isnumeric(x)&&x>=0));
p.addParameter('MinModulation',0,@(x)isscalar(x)&&isnumeric(x)&&x>=0);
p.parse(varargin{:});
opt = p.Results;
gridtype = lower(char(opt.Grid));
if ~ismember(gridtype,{'half','full'})
    error('%s: ''Grid'' must be ''half'' or ''full'', not ''%s''.',mfilename,gridtype);
end

hos = hos(1);
order = double(hos.order);
if ~ismember(order,[3 4])
    error('%s: the zero-lag sum is defined for orders 3 and 4; hos.order is %d.',mfilename,order);
end
if hos.diagonal_slice
    warning('zerolag_integrator:diagonalSlice',...
        ['%s: hos.diagonal_slice is set, so only the point nu = f2-f1 of each line was ' ...
         'estimated (one coefficient per cell, |f1-f2| <= slowpass). Refit with diagonal_slice = false.'],mfilename);
end

fx = hos.freqindx;
len = size(fx.Is,1);                       % reduced coefficients (the trailing NaN slot excluded)
remap = fx.remap;
pdc = fx.PDconj;
fs = hos.sampling_rate;
N = hos.fftN;
df = fs/N;
ndim = order-1;                            % square-form dimensions
lowpass = hos.lowpass;
if isscalar(lowpass)
    lowpass = lowpass*ones(1,ndim);
end
tol = 1e-6;

%%% Square-form axes -> position lookups. Bfreqs{k} are the DFT frequencies
%%% within the lowpass of axis k, in fft order; remap is indexed by them.
pos = cell(1,ndim);
for k = 1:ndim
    kb = round(fx.Bfreqs{k}(:)/df);
    if numel(kb) ~= size(remap,k)
        error('%s: freqindx.Bfreqs{%d} has %d entries but freqindx.remap has %d along axis %d.',...
            mfilename,k,numel(kb),size(remap,k),k);
    end
    pos{k} = zeros(N,1);
    pos{k}(mod(kb,N)+1) = 1:numel(kb);
end

%%% Modulation frequencies, in bins
numax = opt.MaxModulation;
if isempty(numax)
    switch order
        case 4
            % pair sums are kept up to slowpass - passtol (find_principal_domain)
            nmax = floor((min(hos.slowpass) - 0.25*fs/hos.buffersize)/df + tol);
        case 3
            xlp = min(hos.xlowpass);
            if xlp < min(lowpass)
                nmax = ceil(xlp/df - tol) - 1;          % |w3| < xlowpass, strictly
            else
                nmax = floor(min(lowpass)/df + tol);
            end
    end
else
    nmax = floor(numax/df + tol);
end
nmin = ceil(opt.MinModulation/df - tol);
nu = unique([-nmax:-nmin, nmin:nmax]);
nnu = numel(nu);

%%% Carriers on the half grid: f = a*df/2, a = 0..amax
amax = floor(2*min(lowpass)/df + tol);
na = amax+1;
fchalf = (0:amax)*df/2;

%%% Enumerate the (carrier[s], nu) lattice; bins in units of df
switch order
    case 4
        [A,Bc,Nn] = ndgrid(0:amax,0:amax,nu);
        A = A(:); Bc = Bc(:); Nn = Nn(:);
        ok = mod(A+Nn,2)==0 & mod(Bc+Nn,2)==0;    % same lattice: the bins are integers
        A = A(ok); Bc = Bc(ok); Nn = Nn(ok);
        K = [(A+Nn)/2, (Nn-A)/2, -(Bc+Nn)/2, (Bc-Nn)/2];   % f1+nu/2, -f1+nu/2, -f2-nu/2, f2-nu/2
        nrows = na^2;
        row = sub2ind([na na],A+1,Bc+1);
    case 3
        [A,Nn] = ndgrid(0:amax,nu);
        A = A(:); Nn = Nn(:);
        ok = mod(A+Nn,2)==0;
        A = A(ok); Nn = Nn(ok);
        K = [(A+Nn)/2, (Nn-A)/2, -Nn];                       % f+nu/2, -f+nu/2, -nu
        nrows = na;
        row = A+1;
end
kmin = -floor(N/2); kmax = ceil(N/2)-1;    % bins fftfreq(N) represents (guards aliasing)
ok = all(K >= kmin & K <= kmax,2);
row = row(ok); Nn = Nn(ok); K = K(ok,:);
subs = cell(1,ndim);
ok = true(size(row));
for k = 1:ndim
    subs{k} = pos{k}(mod(K(:,k),N)+1);
    ok = ok & subs{k}>0;                   % inside the square form
end
subs = cellfun(@(s)s(ok),subs,'uniformoutput',false);
lin = sub2ind(size(remap),subs{:});
j = double(remap(lin));
c = pdc(lin);
row = row(ok); Nn = Nn(ok);
keep = j <= len;                           % remap == len+1 marks the masked / absent coefficients
row = row(keep); j = j(keep); c = c(keep); Nn = Nn(keep);

I = sparse(row(~c),j(~c),1,nrows,len);
Iconj = sparse(row(c),j(c),1,nrows,len);
[~,nuidx] = ismember(Nn,nu);
present = sparse(row,nuidx,1,nrows,nnu);
presenthalf = present;
fc = fchalf;
n = na;
F = [];

if strcmp(gridtype,'full')
    %%% Linear (order 3) / bilinear (order 4) fold of the half grid onto the DFT grid
    nf = floor(min(lowpass)/df + tol)+1;
    a = (0:amax)';
    switch order
        case 4
            [a,b] = ndgrid(0:amax,0:amax);
            a = a(:); b = b(:);
            ev = mod(a,2)==0 & mod(b,2)==0;
            od = mod(a,2)==1 & mod(b,2)==1;
            r = sub2ind([nf nf],a(ev)/2+1,b(ev)/2+1);
            cidx = sub2ind([na na],a(ev)+1,b(ev)+1);
            w = ones(size(r));
            ao = a(od); bo = b(od);
            co = sub2ind([na na],ao+1,bo+1);
            for da = [-1 1]
                for db = [-1 1]
                    ia = (ao+da)/2; ib = (bo+db)/2;
                    valid = ia>=0 & ia<=nf-1 & ib>=0 & ib<=nf-1;
                    r = [r; sub2ind([nf nf],ia(valid)+1,ib(valid)+1)]; %#ok<AGROW>
                    cidx = [cidx; co(valid)]; %#ok<AGROW>
                    w = [w; 0.25*ones(nnz(valid),1)]; %#ok<AGROW>
                end
            end
            F = sparse(r,cidx,w,nf^2,na^2);
        case 3
            ev = mod(a,2)==0;
            od = ~ev;
            r = a(ev)/2+1;
            cidx = a(ev)+1;
            w = ones(size(r));
            ao = a(od);
            for da = [-1 1]
                ia = (ao+da)/2;
                valid = ia>=0 & ia<=nf-1;
                r = [r; ia(valid)+1]; %#ok<AGROW>
                cidx = [cidx; ao(valid)+1]; %#ok<AGROW>
                w = [w; 0.5*ones(nnz(valid),1)]; %#ok<AGROW>
            end
            F = sparse(r,cidx,w,nf,na);
    end
    I = F*I;
    Iconj = F*Iconj;
    presenthalf = present;
    present = F*present;
    fc = (0:nf-1)*df;
    n = nf;
end

count = full(sum(I,2)+sum(Iconj,2));
if order == 4
    count = reshape(count,n,n);
end
info.fc = fc;
info.n = n;
info.count = count;
info.nu = nu*df;
info.nupresent = unique(Nn)'*df;
info.present = present;
info.presenthalf = presenthalf;
info.fold = F;
info.fchalf = fchalf;
info.grid = gridtype;
info.order = order;
info.len = len;
info.tube = numel(unique(j));
info.df = df;
