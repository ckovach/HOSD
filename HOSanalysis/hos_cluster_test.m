
function stat = hos_cluster_test(hos,yin,xin,nperm,pthresh)

% stat = mdg_cluster_test(hos,yin,[xin],[nperm],[pthresh])
% Implements a cluster-based permutation test for the regression model
% that relates HOS of yin to xin. The default model is intercept-only,
% (equivalent to xin = []) and is thus a test on whether the mean HOS coefficient deviates
% from 0 (in the complex plane). When xin is non-constant, the regression
% model relates the HOS of yin to the weighted average of xin over analysis
% windows (weighted by the analysis window taper).
%
% Inputs:
%   yin: time series whos HOS is the dependent measure
%   xin: time series that gives the independent measure (default is
%        intercept only)
%   nperm: number of permutations (default is 1000). 
%   pthresh: threshold on the p-value from the parametric test, used for
%           cluster detection (default is 0.05. 

if nargin<3
    xin = [];
end

if nargin < 4 || isempty(nperm)
    nperm = 1000;
end
if nargin < 5 || isempty(pthresh)
    pthresh = .05;
end

n = length(hos.B);
if hos.diagonal_slice
    mgi = round(hos.modgram((1:n)'));
    cm = connectmat(true(size(mgi)));
else
    indx = 1:n;
    mgi = round(indx(hos.fullmap));
    mgi = mgi(hos.freqindx.principal_domain>0);
    cm = connectmat(hos.freqindx.principal_domain>0);
end


Adj = spalloc(n,n,n*8);
Adj(mgi(:),mgi(:)) = cm.Adj;

[out,FFY,x,reg_args] = hos.hos_regress(yin,xin);

thri = out.pval<pthresh;
pAdj = Adj(thri,thri);
deg = zeros(size(thri));
deg(thri) = sum(pAdj);
thri = thri & deg>0;

gr = graph(Adj(thri,thri));

stat = out;

stat.clusters = zeros(n,1);
[stat.clusters(thri ),stat.clsize] = conncomp(gr);
for k = 1:length(stat.clsize)
    stat.sdev(k) = sum(out.dev(stat.clusters==k));
    stat.maxdev(k) = max(out.dev(stat.clusters==k));
    stat.minp(k) = min(out.pval(stat.clusters==k));
end

%%
fprintf('\nPermutation: ',0);
nfp=0;
for pidx = 1:nperm+1
      %Permute by randomizing phase across windows. 
      if pidx == 1
          phrand = 1;
      else
          phrand = exp(1i*2*pi*rand(size(FFY,1),1));
      end
      [outp.beta,outp.dev,outp.pval] = complexglm(FFY.*phrand,x,'diagonly',false,reg_args{:});
      statp(pidx).clusters = zeros(n,1);
        

      thri = outp.pval<pthresh;
    pAdj = Adj(thri,thri);
    deg = zeros(size(thri));
    deg(thri) = sum(pAdj);
    thri = thri & deg>0;

        
        if sum(thri)==0
            statp(pidx)= struct('clusters',[],...
                 'sdev',[],...
                 'maxdev',[],...
                 'minp',[],...
                 'maxclsize',double(any(outp.pval<pthresh)),...
                 'maxsdev', max(outp.dev),...
                 'maxmaxdev',max(outp.dev),...
                 'clsize',[]);
        else
            
            gr = graph(Adj(thri,thri));
    
            [statp(pidx).clusters(thri ),statp(pidx).clsize] = conncomp(gr);
            statp(pidx).maxclsize = max(statp(pidx).clsize);
            for k = 1:length(statp(pidx).clsize)
                statp(pidx).sdev(k) = sum(outp.dev(statp(pidx).clusters==k));
                statp(pidx).maxdev(k) = max(outp.dev(statp(pidx).clusters==k));
                statp(pidx).minp(k) = min(outp.pval(statp(pidx).clusters==k));
            end
             statp(pidx).maxsdev = max( statp(pidx).sdev);
             statp(pidx).maxmaxdev = max( statp(pidx).maxdev);
        end
     nfp = fprintf([repmat('\b',1,nfp),'%i'],pidx)-nfp;
end

stat.cluster_sdev_pval = min(1 - mean(statp(1).sdev>[statp(2:end).maxsdev]') + 1./nperm/2,1);
stat.cluster_size_pval = min(1 - mean(statp(1).clsize>[statp(2:end).maxclsize]') + 1./nperm/2,1);
stat.cluster_maxdev_pval = min(1 - mean(statp(1).maxdev>[statp(2:end).maxmaxdev]') + 1./nperm/2,1);

