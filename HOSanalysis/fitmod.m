function out = fitmod(mdl)


% Fit the model

% do_llr_tests = false;
use_glmfit = false;
out.model = mdl;

regs = mdl.designMtx;


X =[regs.value];

if use_glmfit
    [b,devfull,stat] = glmfit(X,mdl.response,mdl.modelType); 
	out.b = b([2:end,1]);
else
    X(:,end+1)=1;
    [b,H,LL] = vectorglm(X,mdl.response,[],mdl.modelType,'gaussreg',1e-6); 
    devfull = -2*LL;
    stat.covb = -H^-1;
    stat.beta = b;
    out.b = b;
    
end
codes = [regs.code];

out.devfull = devfull;
out.stat = stat;
if use_glmfit
    out.yfit = [ones(size(X,1),1),X]*b;
     out.ysd = sum(X.*(X*stat.covb(2:end,2:end)));
    out.intercept = b(1);
else
    out.yfit = X*b;
    out.ysd = sum(X.*(X*stat.covb));
    out.intercept = b(end);
end
out.aic = devfull + 2*length(b);
out.bic = devfull + length(b)*log(size(X,1));

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
    regs(regi).waldpval = 1-chi2cdf(full(waldstat),length(bsub));
    
    if ismember(k,[mdl.do_llr_tests{single_reg_tests}])
        
        if use_glmfit
             [bexcl,devred] = glmfit(X(: ,find([regs.codevec]~=codes(k))),mdl.response,mdl.modelType); 
             bexcl=bexcl([2:end 1]);
        else
            [bexcl,~,LL] = vectorglm(X(: ,[find([regs.codevec]~=codes(k)),end]),mdl.response,[],mdl.modelType);
            devred = full(-2*LL);
        end    
       
        regs(regi).llrpval= 1-chi2cdf(devred-devfull,length(subi));
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
            [bexcl,~,LL] = vectorglm(X(: ,[getreg,end]),mdl.response,[],mdl.modelType);
            devred = full(-2*LL);
        end    
        
        subi=find(ismember([regs.codevec],codes(mdl.do_llr_tests{kk})));
        llrpval= 1-chi2cdf(devred-devfull,length(subi));
        ddev= devred-devfull;
%         bexcl=bexcl;
        out.llrtests(kk) = struct('llrpval', llrpval,'ddev', ddev,'bexcl', bexcl,'regs',mdl.do_llr_tests{kk});
%         out.llrtestss(kk).regs = getreg;
  
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
