function out = fitmod(mdl)


% Fit the model

do_llr_tests = false;
use_glmfit = false;
out.model = mdl;

regs = mdl.designMtx;


X =[regs.value];

if use_glmfit
    [b,devfull,stat] = glmfit(X,mdl.response,mdl.modelType); 
else
    X(:,end+1)=1;
    [b,H,LL] = vectorglm(X,mdl.response,[],mdl.modelType,'gaussreg',1e-6); 
    devfull = -2*LL;
    stat.covb = -H^-1;
end
codes = [regs.code];

out.intercept = b(1);
out.b = b;
out.devfull = devfull;
out.stat = stat;
if use_glmfit
    out.yfit = [ones(size(X,1),1),X]*b;
     out.ysd = sum(X.*(X*stat.covb(2:end,2:end)));
else
    out.yfit = X*b;
    out.ysd = sum(X.*(X*stat.covb));
end
out.aic = devfull + 2*length(b);
out.bic = devfull + length(b)*log(size(X,1));

for k = 1:length(codes)

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
    regs(regi).waldpval = 1-chi2cdf(full(waldstat),length(bsub));
        
    if do_llr_tests
        
        if use_glmfit
             [~,devred] = glmfit(X(: ,[regs.codevec]~=codes(k)),mdl.response,mdl.modelType); 
        else
            [~,~,LL] = vectorglm(X(: ,[regs.codevec]~=codes(k)),mdl.response,[],mdl.modelType);
            devred = full(-2*LL);
        end    
        regs(regi).llrpval= 1-chi2cdf(devred-devfull,length(subi));
        regs(regi).ddev= devred-devfull;
    end

end
switch mdl.modelType
    case 'binomial'
        pfun=@(x)1./(1+exp(-x));
end

csy = cumsum(pfun(out.yfit))/sum(mdl.response);
try
    [~,out.kstest] = kstest(find(mdl.response),[find(mdl.response),csy(find(mdl.response))]);
catch
    out.kstest = nan;
end
for k = 1:length(regs)
    regs(k).value(:) = [];
end
out.regressors = regs;
