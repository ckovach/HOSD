 function out = hos_regress(me,yin,xin,do_permtest,varargin)

% out = hos_regress(me,Y,[X],[do_permtest])
% Linear regression on each element of the HOS array
% computed for Y with regressors in X using a parametric complex Gaussian model.
% 
% Inputs:
%    Y - dependent variable on which HOS statistics are to be computed. Y
%    may be either a column vector or me.buffersize x N matrix. 
%    X - Independent variable. If not given, X is unit valued, so that the
%        test is on the mean. X may be a size(Y,1) x P matrix or a N x P
%        matrix. If the former, then X is converted to an N x P matrix by
%        taking the weighted average value of each chopped segment, where the
%        windowing specified in me.window provides the weighting.
%    do_permtest - carry out a permutation test by shuffling the values of X
%        and generating an empirical distribution for deviance under
%        shuffling. Default is false.
%
% Outputs:
%   out - A structure with fields:
%        .beta - Beta (P+1) x M array of fitted beta values and intercept in the regression model
%     	 .dev  - 1 x M array of deviance.
%        .pval - 1 x M array of P-values for the parametric complex Gaussian model. 
%        .sigma - 1 x M array of error variance estimates
%        .iXX - inverse covariance of X (error covariance for .beta(:,k) is
%               iXX*sigma(k).
%        .wald - (P+1) x M array of Wald statistics for each parameter 
%                estimate (wald(:,k) = beta(:,k)./sqrt(diag(iXX)*sigma(k)) )
%        .permp - 1 x M array of permutation p-values if do_permtest is
%        true.
%
% see COMPLEXGLM HOSOBJECT/SETREG
%
% Copyright Christopher K. Kovach, University of Iowa 2018-2021

    if nargin < 4
        do_permtest = false; %Do a permutation test to verify significant results
    end
    if  nargin > 3 && do_permtest > 1
        maxpermn = dopermtest;
    else
        maxpermn = 5e5;
    end
    
    reg_args = {};
    if nargin < 3 || isempty(xin)
        if size(yin,2) == 1
           xin = ones(size(yin,1),1);
        else
            xin = ones(size(yin,2),1);
        end
        reg_args = [reg_args,{'intercept',false}];
    end
    me(1).regressor = xin;

    if size(yin,2)==1
        Ychop = me(1).chop_input(yin,true);
        keep = ~any(isnan(Ychop)); %Discard any samples containing nans
        Ychop = Ychop(:,keep,:);
        if size(xin,1)==size(yin,1)
            for k = 1:size(xin,2)
                Xchop = me(1).chop_input(xin(:,k),true);
                Xchop = Xchop(:,keep,:);
                x(:,k) = sum(Xchop)/sum(me(1).win);                
            end
        else
            x = xin(keep,:);
        end
        FY = fft(Ychop);
    else
        FY = fft(yin);
        x = xin;
        
    end
    FFY = conj(FY(me(1).freqindx.Is(:,me(1).order),:));          
    for k = me(1).order-1:-1:1
        FYk = FY(me(1).freqindx.Is(:,k),:);
        FFY = FFY.*FYk;
    end
%             [~,out.dev0] = complexglm(FFX',ones(size(x,1),1),'diagonly',false,'intercept',false);
    [out.beta,out.dev,out.pval,out.iXX,out.sigma] = complexglm(FFY',x,'diagonly',false,reg_args{:});
    out.beta(:,end+1) = nan;
    out.dev(end+1) = nan;
    out.pval(end+1)=nan;

    out.sigma(end+1)=nan;
    me(1).regstat=out;
    se = sqrt(diag(out.iXX)*out.sigma);
    out.wald = out.beta./se;
    if do_permtest
       permi = 1;
       den = zeros(size(out.pval(1:end-1)));
       num = den;
       num_threshold = [5 4 3 2 1]; % Number of super-threshold permutations at which to stop at ceil(log10(permutation_index)
       fp = fprintf('\nPermutation %i, Nsig: %i',0,0);
       %%% Permutation test with number of permutations adjusted
       %%% to nominal p-value. 
       reseed;
       geti=1;
       while any(out.pval<=1/(permi*2)) && permi<maxpermn && sum(geti)>0
           rp = randperm(size(FFY,2));
           if mod(permi,10.^floor(log10(permi)-2))==0
           geti = num < num_threshold(min(end,ceil(log10(permi+1))));
               fp = fprintf([repmat('\b',1,fp),'Permutation %i, Nsig: %i'],permi,sum(geti))-fp;
           end
           [~,dev] = complexglm(FFY(geti,rp)',x,'diagonly',false,reg_args{:});
           den(geti) = den(geti)+1;
           num(geti) =  num(geti)+(out.dev(geti)<dev);
           permi=permi+1;
       end
       out.permp = num./den;
       out.permp(end+1) = nan;
%                out.nperms = den;

    end

end