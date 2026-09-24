function S = zerolag_pairmean(FX,hos,info)

% S = zerolag_pairmean(FX,hos,info)
%
% Second-order terms for ZEROLAG_BANDCORR: the mean pair products
%
%    S.m(f,nu) = E[M_f(nu)] = mean over windows of X(f+nu/2) X*(f-nu/2)
%
% on the half-grid carriers info.fchalf and the modulation frequencies
% info.nu of ZEROLAG_INTEGRATOR (NaN where f and nu are on different
% lattices), and S.xbar(nu) = mean over windows of X(nu), from the windowed
% spectra FX (fftN x windows, the matrix handed to
% hosobject/update_bispectrum). The lag phase shift update_bispectrum
% applies (hos.lag) is applied here too, and the windows are weighted
% uniformly, as a batch fit (update_bispectrum(FX,true)) weights them.
%
% Why: with a data window w(t) the mean WVD row of even a stationary signal
% follows w(t)^2, so every pair of rows shares that profile and the order-4
% moment sum reports it as co-modulation (about 0.3 between independent
% bands under a sasaki window). sum_nu E[M_f1(nu)] E[M_f2(nu)]* is that
% term -- the second-order part of the fourth-order cumulant; the remaining
% parts vanish for well-separated positive carriers -- and subtracting it
% leaves the zero-lag COVARIANCE of the rows, so that 'corr' is a genuine
% correlation coefficient. At order 3 the analogous term is
% E[M_f(nu)] E[X(nu)]*, zero for a zero-mean stationary signal.
%
% See also ZEROLAG_INTEGRATOR, ZEROLAG_BANDCORR, HOSOBJECT/UPDATE_BISPECTRUM.

% C. Kovach 2026

hos = hos(1);
N = hos.fftN;
fs = hos.sampling_rate;
df = fs/N;
if size(FX,1) ~= N
    error('%s: FX must have fftN = %d rows (the windowed spectra).',mfilename,N);
end
%%% The lag adjustment update_bispectrum applies
dt = atan2(imag(hos.lag),real(hos.lag))/(2*pi)*N;
delt = hos.radw(:)*dt;
delt(isnan(delt)) = 0;
FX = exp(-1i*delt).*FX;

a = round(2*info.fchalf(:)/df);
nn = round(info.nu(:)'/df);
[A,Nn] = ndgrid(a,nn);
k1 = (A+Nn)/2;
k2 = (A-Nn)/2;
kmin = -floor(N/2); kmax = ceil(N/2)-1;
valid = mod(A+Nn,2)==0 & k1>=kmin & k1<=kmax & k2>=kmin & k2<=kmax;
m = nan(size(A));
b1 = mod(k1(valid),N)+1;
b2 = mod(k2(valid),N)+1;
m(valid) = mean(FX(b1,:).*conj(FX(b2,:)),2);
S.m = m;
S.xbar = mean(FX(mod(nn,N)+1,:),2).';
S.fchalf = info.fchalf;
S.nu = info.nu;
