function [out,X] = fitmod(mdl)


% Fit the model
%
%   out = fitmod(mdl)
%   [out,X] = fitmod(mdl)   also returns the assembled design matrix, with
%                           the constant column appended as its last column.
%
%   mdl.observation_weight (see model.observationWeights) weights each
%   observation's contribution to the log-likelihood. A weight of 0
%   WITHHOLDS that observation: every column of its design row is zeroed,
%   the appended constant column included, and the row is dropped from the
%   likelihood, so it contributes nothing to the score, to the Hessian or
%   to the deviance. That is exact row deletion, but the row stays in place
%   so yfit / ysd / yhat keep one entry per observation. Because the
%   intercept is zero there too, a withheld bin is not read as "a bin at
%   the reference level" and is not pooled into the intercept or into the
%   reference level of a categorical block.
%
%   The fit statistics (devfull, aic, bic, kstest) and the nested LLR
%   refits are computed over the INCLUDED observations only, and bic counts
%   those rather than all rows. out.included, out.nobs and out.nexcluded
%   report what was used.
%
%   Caveat: on a withheld row yfit and ysd are 0 and yhat is pfun(0) --
%   1 for a Poisson model, 0.5 for a binomial one. Those are the arithmetic
%   of an all-zero design row, not predictions; select with out.included
%   before summing or plotting them.

% do_llr_tests = false;
use_glmfit = false;
out.model = mdl;

glmtype = mdl.modelType;
if strcmpi(glmtype,'count')
    glmtype = 'poisson';
end

regs = mdl.designMtx;


X =[regs.value];

w = mdl.observationWeights(); % [] when every observation counts equally
if isempty(w)
    included = true(size(X,1),1);
    glmargs = {};
else
    included = w > 0;
    glmargs = {'weights',w};
end
nobs = sum(included); % observations that enter the estimating equations

if use_glmfit
    if ~isempty(w)
        error('fitmod:observationWeight',...
              'observation_weight needs the vectorglm path: glmfit appends its own intercept, which cannot be zeroed on a withheld row.');
    end
    [b,devfull,stat] = glmfit(X,mdl.response,mdl.modelType);
	out.b = b([2:end,1]);
else
    X(:,end+1)=1;
    if ~isempty(w)
        X(~included,:) = 0; % a withheld observation is zero in EVERY column, the appended constant included
    end
    [b,H,LL,msg] = vectorglm(X,mdl.response,[],glmtype,'gaussreg',1e-6,'showiter',false,glmargs{:});
    devfull = -2*LL;
    stat.covb = -H^-1;
    stat.beta = b;
    out.b = b;
    out.msg = msg; % convergence / conditioning flags from vectorglm
    
end
codes = [regs.code];

out.devfull = devfull;
out.stat = stat;
if use_glmfit
    out.yfit = [ones(size(X,1),1),X]*b;
     out.ysd = sqrt(sum(X.*(X*stat.covb(2:end,2:end)),2)); % std. err. of the linear predictor per sample
    out.intercept = b(1);
else
    out.yfit = X*b;
    out.ysd = sqrt(sum(X.*(X*stat.covb),2)); % std. err. of the linear predictor per sample
    out.intercept = b(end);
end
out.aic = devfull + 2*length(b);
out.bic = devfull + length(b)*log(nobs); % nobs = the observations that were actually fit
out.included = included;                 % logical, one entry per observation
out.nobs = nobs;
out.nexcluded = numel(included) - nobs;  % observations withheld by observation_weight
out.observation_weight = w;              % [] when nothing was withheld or reweighted

if islogical(mdl.do_llr_tests) 
    if mdl.do_llr_tests
        mdl.do_llr_tests = num2cell(1:length(regs));        
    else
        mdl.do_llr_tests = {};
    end
elseif isnumeric(mdl.do_llr_tests)
    mdl.do_llr_tests = num2cell(mdl.do_llr_tests);
end
single_reg_tests = cellfun(@(x)length(x)==1,mdl.do_llr_tests);
for k = 1:length(codes)
%%
    if use_glmfit  
        subi=find([regs.codevec]==codes(k))+1;
    else
        subi=find([regs.codevec]==codes(k));
    end
    bsub = b(subi);
    covbsub = stat.covb(subi,subi);
    regi = find(codes==codes(k));
    regs(regi).beta = bsub;
    regs(regi).covb = covbsub;
    
    if ~isempty(regs(regi).window)
        evw = mdl.get_event_window(regs(regi).window);
        
        [unqlev,~,unqlevi] = unique(regs(regi).levmat(1:end-1,:)','rows');
        regs(regi).windowest=struct('intensity',[],'sd',[],'wald',[],'tt',[]);
        for kk = 1:size(unqlev,1)
            levi = unqlevi==kk;
            regs(regi).windowest.intensity(:,kk) = evw.P*bsub(levi); 
            regs(regi).windowest.sd(:,kk) = sqrt(sum(evw.P.*(evw.P*covbsub(levi,levi)),2)); 
            regs(regi).windowest.wald(:,kk) = regs(regi).windowest.intensity(:,kk)./regs(regi).windowest.sd(:,kk);
        end
        regs(regi).windowest.tt=evw.tt;
    end
    waldstat = bsub'*covbsub^-1*bsub;
    regs(regi).waldstat = waldstat;
    regs(regi).waldpval = gammainc(full(waldstat)/2,length(bsub)/2,'upper'); % chi2 upper tail: no stats toolbox, exact for small p
    
    if ismember(k,[mdl.do_llr_tests{single_reg_tests}])
        
        if use_glmfit
             [bexcl,devred] = glmfit(X(: ,find([regs.codevec]~=codes(k))),mdl.response,mdl.modelType); 
             bexcl=bexcl([2:end 1]);
        else
            [bexcl,~,LL] = vectorglm(X(: ,[find([regs.codevec]~=codes(k)),end]),mdl.response,[],glmtype,'gaussreg',1e-6,'showiter',false,glmargs{:}); % same penalty and the same observations as the full fit so the models nest
            devred = full(-2*LL);
        end    
       
        regs(regi).llrpval= gammainc(max(devred-devfull,0)/2,length(subi)/2,'upper');
        regs(regi).ddev= devred-devfull;
        regs(regi).bexcl=bexcl;
        geti = cellfun(@(x)isequal(k,x),mdl.do_llr_tests);
        out.llrtests(geti) = struct('llrpval', regs(regi).llrpval,'ddev', regs(regi).ddev,'bexcl', regs(regi).bexcl,'regs', mdl.do_llr_tests{geti});
    end

end

for kk = find(~single_reg_tests)   
      getreg = find(~ismember([regs.codevec],codes(mdl.do_llr_tests{kk})));
       if use_glmfit
             [bexcl,devred] = glmfit(X(: ,[1,getreg]),mdl.response,mdl.modelType); 
             bexcl=bexcl([2:end 1]);
        else
            [bexcl,~,LL] = vectorglm(X(: ,[getreg,end]),mdl.response,[],glmtype,'gaussreg',1e-6,'showiter',false,glmargs{:}); % same penalty and the same observations as the full fit so the models nest
            devred = full(-2*LL);
        end    
        
        subi=find(ismember([regs.codevec],codes(mdl.do_llr_tests{kk})));
        llrpval= gammainc(max(devred-devfull,0)/2,length(subi)/2,'upper');
        ddev= devred-devfull;
%         bexcl=bexcl;
        out.llrtests(kk) = struct('llrpval', llrpval,'ddev', ddev,'bexcl', bexcl,'regs',mdl.do_llr_tests{kk});
%         out.llrtestss(kk).regs = getreg;
  
end

switch lower(glmtype)
    case 'binomial'
        pfun=@(x)1./(1+exp(-x));
    case 'poisson'
        pfun=@(x)exp(x);
    otherwise
        error('Unrecognized model type %s',mdl.modelType);
end
out.yhat = pfun(out.yfit); % expected count (poisson) / event probability (binomial) per sample

% Time-rescaling KS test over the included observations only: a withheld row
% has no fitted intensity, so it must not accumulate one.
yhinc  = out.yhat(included);
respinc = mdl.response(included);
csy = cumsum(yhinc)/sum(respinc);
evi = find(respinc);
try
    [~,out.kstest] = kstest(evi,[evi,csy(evi)]);
catch
    out.kstest = nan;
end
for k = 1:length(regs)
    regs(k).value(:) = [];
end
out.regressors = regs;
