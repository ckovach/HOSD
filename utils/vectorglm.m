
function [theta,H,LL,msg] = vectorglm(X,y,ntrials,varargin)

% Fast vectorized fitting of multiple poisson or binomial logistic 
% regressions.
% 
% 
% [theta,H,LL,msg] = vectorglm(X,Y,ntrials)
% 
%  X : N x K design matrix
%  
%  Y: a matrix of event counts. Note: for binomial models, if any Y > 1,
%     then you must inclued the ...,'Ntot',N,... argument containing the total
%     number of trials for each row of X. Y(i) must be <= N(i).
%
%  nrows : 1 X M vector where ntrials(i) gives the number of contiguous rows in X belonging to
%           the ith regression.
%
%  theta : MLE parameter estimates
%  H  :  Hessian matrix
%  LL :  Log-likelihood.
%  msg: struct with fields
%         .badcond: true if hessian is badly conditioned or singular
%         .converged: false if estimate did not converge before the maximum
%                   number of iterations
%         .Niter: = number of iterations.
%
%  [...] = vectorlogistic(X,ntrials,options)
%
%    options:
%       'binomial' : fit binomial logistic model (default)
%       'Poisson' : fit poisson model (default)
%       'Laplreg', L : %Specify laplacian regularization, sqrt(l)*|r| (the
%                    SQUARE of the weighting on the norm), L can be scalar
%                    vector or matrix.
%       'lasso', L  :  lasso regularization
%       'gaussreg',H : Gaussian regularization. H is the 'precision' matrix
%                   for the gaussian prior (inverse of
%                   variance-covariance).
%       'weights',W : Observation weights, one non-negative value per row of
%                   Y. The log-likelihood, score and Hessian are summed as
%                     LL  =  sum_i W_i * l_i(theta)
%                     DLL =  sum_i W_i * x_i * (y_i - mu_i)
%                     D2LL= -sum_i W_i * v_i * x_i * x_i'
%                   so W_i = 0 withholds observation i (it has no influence
%                   on model fit). An integer W_i is equivalent to repeating
%                   that observation W_i times.Default: every
%                   observation weighted 1.
%
%     For other options, see the menu in the m-file

%
% C Kovach 2011
%


type = 'binomial';



% L1reg = 0;
% L2reg = 0;
Hreg = sparse(0);
Lreg = sparse(0);
L1reg = sparse(0);
fix  = 0;
diagHess = false;
Ntot = ones(size(y));
obsweight = []; % observation weights ([] = every observation weighted 1)
maxiter = 200;
maxlasso = 10; % Maximum number of times to repeat lasso procedure.
showiter = true;
tol = 1e-6;
printiter = 10; %print information every i'th iteration.
LC = [];
do_grcheck = false;
initTheta = [];
i = 1;
runiter = true;

%%%% Options are described below %%%%
while i <= length(varargin)
    
    switch lower(varargin{i})
        case {'poisson'}
            type = 'poisson';
            
        case {'binomial'}
            type = 'binomial';
         case {'gaussreg','regularization','l2reg'}   %Specify gaussian regularization  the SQUARE of the weighting on the L2 norm (sqrt(g)*|r|^2)
            Hreg = varargin{i+1};
            i = i+1;            
%        case {'l2reg', 'ridge'}   %Specify ridge regularization (L2, ridge), this is like gaussreg but with weighting that increases in proprtion to the data.
%             L2reg = varargin{i+1};
%             i = i+1;                 
        case {'laplreg'}   %Specify laplacian regularization, sqrt(l)*|r| (the SQUARE of the weighting on the norm
            Lreg = varargin{i+1};
            i = i+1;      
        case {'l1reg','lasso'}   %Specify laplacian regularization, sqrt(l)*|r| (the SQUARE of the weighting on the norm
            L1reg = varargin{i+1};
            i = i+1;      
%         case {'lasso','l1reg'}   %Lasso estimator, which is proportional to the quantity of data.
%             L1reg = varargin{i+1};
%             i = i+1; 
        case {'inittheta'}
            initTheta = varargin{i+1};
            i = i+1;      
       case {'ntot'}
            Ntot = varargin{i+1};
            i = i+1;
       case {'weights','weight','obsweight','observation_weight'}   %Non-negative weight on each observation's log-likelihood; 0 withholds the observation
            obsweight = varargin{i+1};
            i = i+1;
       case {'diag'}
            diagHess = varargin{i+1};
            i = i+1;      
       case {'maxiter'}
            maxiter = varargin{i+1};
            i = i+1;      
       case {'tol','tolerance'}
            tol = varargin{i+1};
            i = i+1;      
       case {'runiter'}
            runiter = varargin{i+1};
            i = i+1;      
       case {'grcheck'}
            do_grcheck = varargin{i+1};
            i = i+1;      
       case {'showiter','showprog'}
            showiter = varargin{i+1};
            i = i+1;      
     case 'linearconstraint'     %a matrix of linear constraints on maximization, which is a contrast matrix such that
              LC = varargin{i+1};   %  a*C = 0 if c is the same length as
              i = i+1;              %  the parameters and a*C - C(end,:) = 0
                                    %  if it is longer by 1. Lagrange
                                    %  multipliers are appended to the
                                    % parameter vector.
        case 'fix'
            fix = varargin{i+1};
            i = i+1;  
      
        otherwise
            error('%s is not a recognized keyword',varargin{i})
    end
    i = i+1;
end


if do_grcheck && ~GRcheck
    error('This function needs the gaze reader toolbox to run')
end





if nargin < 3 || isempty(ntrials)
    ntrials = length(y);
end

if length(ntrials) == 1
    ntrials = ones(length(y)./ntrials,1)*ntrials;
end

if isempty(initTheta)
    initTheta = zeros(length(ntrials)*size(X,2),1);
end


spX = sparseblock(X+eps,ntrials,'transpose');

y = sparse(y);

%%% Observation weights. obsw is left as the scalar 1 when no weights were
%%% given, so that every 'obsw.*' below is an exact no-op and the unweighted
%%% path is arithmetically identical to the code before weights existed.
if isempty(obsweight)
    obsw = 1;
else
    if islogical(obsweight)
        obsweight = double(obsweight);
    end
    obsw = full(double(obsweight(:)));
    if isscalar(obsw)
        obsw = obsw*ones(size(y,1),1);
    elseif numel(obsw) ~= size(y,1)
        error('Weights must have one entry per observation (%i), found %i.',size(y,1),numel(obsw))
    end
    if any(~isfinite(obsw)) || any(obsw < 0)
        error('Weights must be finite and non-negative (0 withholds an observation).')
    end
    if all(obsw == 1)
        obsw = 1; % nothing to weight: fall back on the exact no-op
    end
end



if min(size(Hreg)) == 1 %&& ~isequal(Hreg ,0)
    Hreg = diag(sparse(ones(length(initTheta),1).*Hreg(:)));
end


if min(size(Lreg)) == 1 %&& Lreg >0
    Lreg = diag(sparse(ones(length(initTheta),1).*Lreg(:)));
end

Npar = length(initTheta);
%Create constraint matrix for fixed parameters
if any(fix(:))
%     constraint = true;
    fixed = find(fix);
    LC = zeros(Npar,1);
    for i = 1:length(fixed)
        LC(i,fixed(i)) = 1;
    end
        LC(end+1,:) = InitTheta(fixed);
end

if ~isempty(LC)    
    if size(LC,2) == Npar
        LC(:,end+1) = 0;
    end
    if size(LC,2) ~= Npar+1
        error('Contraint matrix must have as many columns as the number of parameters')
    end
%     constraint = 1;
    
    lgm = zeros(size(LC,1),1);
else
    lgm = []; %lagrange multipliers    
%     constraint = 0;
end
nlgm = size(LC,1);


%%% weighted length of a vector
wleng = @(x, M ) sparse(sqrt(sum(x'*M*x)+eps));

%%% A matrix which will be used to compute 0th and 1st  derivatives of
%%% the prior




if any(Lreg(:)~=0)
    RegMat = @(theta,a) ( a*Hreg + Lreg./wleng(theta,Lreg) )*theta ;
    
    laplregfun =@(theta)  Lreg/wleng(theta,Lreg) - ((Lreg*theta)*(Lreg*theta)')/wleng(theta,Lreg)^3 ; %%% Gaussian and laplace prior 

    D2Regfun =@(theta) Hreg + laplregfun(theta); %%% Gaussian prior only
else
%     laplregfun = @(theta) 0;
    RegMat = @(theta,a) a*Hreg*theta ; 
        D2Regfun =@(theta) Hreg ; %%% Gaussian prior only

end    





if all(L1reg == 0)
    l1regfun = @(theta) 0;
else
    if length(L1reg) ==1
        L1reg = ones(size(initTheta))*L1reg;
    end
    L1reg = sparse(L1reg);
    l1regfun =@(theta)  diag(L1reg)*sign(theta); %%% Gaussian and laplace prior 
end

l1regfuno = l1regfun;

if strcmpi(type,'binomial') && any(y>Ntot) % Poisson counts are unbounded
    error('y must be a count which is less than Ntot')
end



theta = initTheta;

recomputeHessian = 20;
% del = Inf;
% discard = false;

levstep = 2;
levmar =  .1;




switch lower(type)
    
    case 'binomial'
        
        rho = @(th) spX*th;
        
        %Probability
        Pfun = @(th) exp( rho(th) )./(1+exp(rho(th)));
        
        %Likelihood function (obsw weights each observation's contribution)
        LLvec = @(th) obsw.*(rho(th).*y - Ntot.*log(1+exp(rho(th))));
        LLfun = @(th) sum(LLvec(th)) - th'*(RegMat(th,1) + L1reg.*sign(th));

        %First derivative of the likelihood function
        DLLfun = @(th)  ((obsw.*(repmat(y,1,size(th,2)) - repmat(Ntot,1,size(th,2)).*Pfun(th)))'*spX)'- RegMat(th,2) - l1regfun(th);

        if ~diagHess
            D2LLfun = @(th)   -spX'*diag( sparse(obsw.*Ntot.*Pfun(th).*(1-Pfun(th))) )*spX - sparse(D2Regfun(th));
        else
            D2LLfun = @(th)    -diag(sparse(obsw.*Ntot.*Pfun(th).*(1-Pfun(th)))'*(spX.*spX)) - sparse(D2Regfun(th));
%             D2LLinv = @(th)    -diag((sparse(Ntot.*Pfun(th).*(1-Pfun(th)))'*(spX.*spX)).^-1) - sparse(D2Regfun(th));
        end

    case 'poisson'
        %INtensity
        rho = @(th) spX*th;
        Pfun = @(th) exp( rho(th) );

        LLvec = @(th) obsw.*(rho(th).*y - Pfun(th));
        LLfun = @(th) sum( LLvec(th) ) - th'*(RegMat(th,1) + L1reg.*sign(th)); % scalar penalised log-likelihood, as for the binomial case
        DLLfun = @(th)  ((obsw.*(y - Pfun(th)))'*spX)'- RegMat(th,2) - l1regfun(th);
        D2LLfun = @(th)   -spX'*diag( sparse(obsw.*Pfun(th)) )*spX - sparse(D2Regfun(th)); % Hessian of the log-likelihood is negative definite
%         DLLfun = @(th) sum((y'*- Pfun(th));
        
        
        
end

DLLfuno = DLLfun;
% LLfuno = LLfun;
D2LLfuno = D2LLfun;
% RegMato = RegMat;

if any(L1reg(:) ~=0)
    discard = DLLfun(theta + 1e-6).*DLLfun(theta - 1e-6) < -1e-6;

    discard = discard & abs(DLLfun(theta.*0)+ l1regfun(theta*0)) < L1reg &...
                sign(DLLfun(theta.*0).*DLLfun(theta))>0;
%   discard = false(size(theta));
%     thetanew(discard) = 0;
    
 
else
    discard = 0;
end

       
thetafull = theta;


spXfull = spX;

discold = nan;

% Lregfull = Lreg;
L1regfull = L1reg;
% ths = [];
ll = [];
nlasso = 0;
while ~isequal(discard,discold) && nlasso < maxlasso

    iter = 0;
    del = Inf;

    
    if ~isempty(LC)
%        DLLfun = @(theta) cat(1,DLLfuno(theta) + LC(1:end-1,:)*lgm,[theta; -1]'*LC); %Derivative of parameters under linear constraint
        DLLfun = @(theta) cat(1,DLLfuno(theta) + LC(:,1:end-1)'*lgm,LC*[theta; -1]); %Derivative of parameters under linear constraint
        D2LLfun =  @(theta) cat(1,...
                    cat(2,D2LLfuno(theta),LC(:,1:end-1)' ) ,...
                    [LC(:,1:end-1),-1e-9*eye(nlgm)]) ;             % Derivative of Lagrange Multipliers;    
    end


    if any(L1regfull(:)~=0)
        
        spX = spXfull(:,~discard);
        theta = thetafull(~discard);
        L1reg = L1regfull(~discard);
       RegMat = @(theta,a)  a*Hreg(~discard,~discard)*theta ;%+ Lreg./wleng(theta,Lreg) )*theta ; 
%         l1regfun =@(theta)  L1reg(~discard).*sign(theta); %%% Gaussian and laplace prior 
        l1regfun =@(th)  diag(L1reg)*sign(th);

        D2Regfun =@(theta) Hreg(~discard,~discard) ; %%% Gaussian prior only

        switch lower(type)

            case 'binomial'

                rho = @(th) spX*th;

                %Probability
                Pfun = @(th) exp( rho(th) )./(1+exp(rho(th)));

                %Likelihood function (obsw weights each observation's contribution)
                LLvec = @(th) obsw.*(rho(th).*y - Ntot.*log(1+exp(rho(th))));
%                 LLfun = @(th) sum(rho(th).*y - Ntot.*log(1+exp(rho(th)))) - th'*(RegMat(th,1) + L1reg(~discard).*sign(th));
                LLfun = @(th) sum(LLvec(th)) - th'*(RegMat(th,1) + L1reg.*sign(th));

                %First derivative of the likelihood function
%                 DLLfun = @(th)  ((y - Ntot.*Pfun(th))'*spX)'- RegMat(th,2) - l1regfun(th);
                DLLfun = @(th)  ((obsw.*(y - Ntot.*Pfun(th)))'*spX)'- RegMat(th,2) - l1regfun(th);

                if ~diagHess
                    D2LLfun = @(th)   -spX'*diag( sparse(obsw.*Ntot.*Pfun(th).*(1-Pfun(th))) )*spX - sparse(D2Regfun(th));
                else
                    D2LLfun = @(th)    -diag(sparse(obsw.*Ntot.*Pfun(th).*(1-Pfun(th)))'*(spX.*spX)) - sparse(D2Regfun(th));
        %             D2LLinv = @(th)    -diag((sparse(Ntot.*Pfun(th).*(1-Pfun(th)))'*(spX.*spX)).^-1) - sparse(D2Regfun(th));
                end

            case 'poisson'
                %INtensity
                rho = @(th) spX*th;
                Pfun = @(th) exp( rho(th) );

                LLvec = @(th) obsw.*(rho(th).*y - Pfun(th)); % keep the per-sample likelihood on the reduced design too
                LLfun = @(th) sum(LLvec(th)) - th'*(RegMat(th,1) + L1reg.*sign(th));
                DLLfun = @(th)  ((obsw.*(y - Pfun(th)))'*spX)'- RegMat(th,2) - l1regfun(th);
                D2LLfun = @(th)   -spX'*diag( sparse(obsw.*Pfun(th)) )*spX - sparse(D2Regfun(th)); % negative definite
        %         DLLfun = @(th) sum((y'*- Pfun(th));



        end     
    end
    

    
%     fixed = find(fix | discard);
%     LC = zeros(Npar,1);
%     for i = 1:length(fixed)
%         LC(i,fixed(i)) = 1;
%     end
%         LC(end+1,:) = InitTheta(fixed);


spey = speye(length(DLLfun(theta)));


    % Newton-raphson gradient ascent
    lastlap = false;
    stop = false;
    
    while  runiter && ~stop 

        if ~diagHess
            if mod(iter,recomputeHessian) == 0 || lastlap
                try  
                       chR = chol(-D2LLfun(theta) + spey*levmar );  %%Using the cholesky factorization is faster than computing the inverse, but might be worse for poorly conditioned matrices
                catch
                    levmar = .01;
                    try
                   chR = chol(-D2LLfun(theta) + spey*levmar );  %%Using the cholesky factorization is faster than computing the inverse, but might be worse for poorly conditioned matrices
                    catch
                        warning('Convergence error')
                        break
                    end
                end
               %                H = (-D2LLfun(theta) + spey*levmar );
%                while condest(H) > 1e5
%                      levmar = max( 1e-3, levmar*2);
%                      H = (-D2LLfun(theta) + spey*levmar );
%                end
            end
           if showiter && mod(iter,printiter) == 0
               fprintf('\n%i  del: %f',iter,del)
           end
%                 dtheta = H\DLLfun(theta); 
            dtheta = chR\(chR'\sparse(DLLfun(theta))); 
    %         chR = D2LLfun(theta);
    %         dtheta = chR^-1*DLLfun(theta);
            
            
    
        else
            dtheta = -( D2LLinv(theta))*sparse(DLLfun(theta));

        end

%     dtheta = dtheta.*(~discard);
    thetanew = theta + dtheta;
%     if L1reg ~=0
%         discard = discard | thetanew.*theta < -1e-6;
%     
%         discard = discard & abs(DLLfun(theta.*0)+ l1regfun(theta*0)) < L1reg &...
%                     sign(DLLfun(theta.*0).*DLLfun(theta))>0;
%     
%         thetanew(discard) = 0;
%         theta(discard) = 0;
%     end
   
          dll = LLfun(thetanew) - LLfun(theta);
          del = full(sum(dtheta.^2))./sum(theta.^2+eps);
%         del = dll;
%         
        if dll >= 0
            theta = thetanew;
            levmar = levmar./levstep;
        else
               levmar = levmar.*levstep + .001;
            try
               chR = chol(-D2LLfun(theta) + spey*levmar);
            catch  %#ok<CTCH>
               levmar = levmar.*levstep + .1;
               chR = chol(-D2LLfun(theta) + spey*levmar);
            end
            
           if showiter
               fprintf('\n%i  del: %f',iter,del)
           end
        end


        iter = iter+1;

        if iter > maxiter
            warning('Maximum iterations reached')
            break
        elseif isnan(del)
            warning('NaN detecting. Stopping.')
            break
        end
        
        if del + sqrt( levmar*sum(del.^2)./sum(theta.^2)) < tol %%% damping dependent term avoids stopping when the damping becomes high
             if lastlap
                 stop = true;
             else
                 lastlap = true;
             end         
         else
             lastlap = false;
        end
         
%         ll(end+1) = LLfun(theta);
%         ths(~discard,end+1) = theta;
    end
    if any( L1regfull(:)~=0)

        discold = discard;
        thetafull(~discard) = theta;
        thetafull(discard) = 0;
    %     fprintf('\n%6i',0);
        thth = repmat(thetafull,1,length(thetafull)) - diag(thetafull);
        dlth = diag(DLLfuno(thth));
        
        %%% IF the sign of the slope changes on either side of 0, then the
        %%% peak is at zero.
        spe =speye(length(thetafull));
        pkval = (dlth+diag(-l1regfuno(thth + 1e-6*spe ))).*...
                  (dlth+diag(-l1regfuno(thth - 1e-6*spe )));
        
        %%% resurrect a feature if pkval exceed some threshold
        discard = pkval < -1e-6 | discold & pkval > .25;
        nlasso = nlasso+1;
    else
        thetafull = theta;
        discard = 0;
        discold = 0;
    end
    
%     if any(L1reg~=0)
%         for i = 1:length(thetafull)
%             zth = thetafull;
%             dzth = 0*thetafull;
%             zth(i) = 0;
%             dzth(i) = 1e-6;
%             dl = DLLfuno(zth + dzth) .* DLLfuno(zth-dzth);
%             discard(i) = dl(i) < 1e-6;
% %             i
%             if mod(i,10) ==0;
%                 fprintf('\b\b\b\b\b\b%6i',i);
%             end
%         end    
%     end
  

end

theta = thetafull;

H = D2LLfuno(theta);

% LL = LLfuno(theta);

cdtol = 1e6;
msg.badcond = condest(H) > cdtol;
msg.converged = iter < maxiter;
msg.Niter = iter;

if ~isequal(discard,0)
    spX = spXfull(:,~discard);
    theta = thetafull(~discard);
    L1reg = L1regfull(~discard);
    RegMat = @(theta,a)  a*Hreg(~discard,~discard)*theta ;%+ Lreg./wleng(theta,Lreg) )*theta ; 
else
    discard = false(size(theta));
end
% LL = sum(sparseblock(LLvec(theta),ntrials,'transpose')) - sum(sparseblock(theta.*(RegMat(theta,1) + L1reg.*sign(theta)),size(spX(:,~discard),2),'transpose'));
RM = 0*thetafull; RM(~discard) = RegMat(theta,1);
LL = sum(sparseblock(LLvec(theta),ntrials,'transpose')) - sum(sparseblock(theta.*(RM(~discard,:) + L1reg.*sign(theta)),size(X(:,~discard),2),'transpose'));








function gr_present = GRcheck

% This function checks if GazeReader (GR) is present in the current path. If
% not, it offers to attempt a fresh installation through svn. It returns
% a value of 1 if GR is present or successfully installed and 0 otherwise.
%
% Subversion can be downloaded at http://subversion.apache.org/packages.html
%

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


%%% URL for current GazeReader repository
grurl = 'https://saccade.neurosurgery.uiowa.edu/GazeReader/0.1';

if ~exist('GazeReader.m','file')
    beep
    
    
   fprintf(['\n--------MISSING FILES---------\nNecessary files appear not to be in the current path.\n\nPlease install GazeReader and add to path\n',...
            '%s\n\n',grurl]) 
    inp = 'x';    
    while ~ismember(lower(inp),{'y','n'})
        inp = input(sprintf(['\nDo you want me to try to install GazeReader now?\n(This requires',...
                    ' command line svn, and username/passwd access)\nY/N:']),'s');
    end
    
    if isequal(lower(inp),'y')
        
        
        fprintf('\nChoose a place to install..')
        installdir = '';
        while exist(fullfile(installdir,'.svn'),'dir') > 0
           installdir = fullfile('..',installdir); 
        end
                
        installdir = uigetdir(installdir,'Choose where to install GazeReader.');
        while exist(fullfile(installdir,'.svn'),'dir') > 0  && ~exist(fullfile(installdir,'GazeReader.m'),'file')
            warndlg('Selected directory must not be a working copy of a subversion repository. Please choose another.')
            installdir = uigetdir(installdir,'Selected directory must not be a subversion repository.');        
            if installdir == 0
                return
            end
        end    
        
        if exist(installdir,'dir') > 0 && exist(fullfile(installdir,'GazeReader.m'),'file')
            fprintf('\nAdding existing version to the matlab path.')
            addpath(installdir)
            gr_present = true;
            return
        elseif exist(installdir,'dir') > 0
            grpath = fullfile(installdir,'GazeReader');
        else
            grpath = installdir;
        end
        
        com = sprintf('svn co %s %s',grurl,grpath);
        
        %%% Attempt to checkout using subversion at the system command-line
        [stat,res] = system(com);  
        
        if stat >0  %%% If command returned an error
            beep
            if ispc                  
                    url = sprintf('download a free version at\n\n\t\thttp://www.sliksvn.com/en/download');
            elseif ismac 
                    url = sprintf('download a free version at\n\n\t\thttp://www.open.collab.net/downloads/community/');
            elseif isunix
                    url = sprintf('install it from your \nrepository (eg run  ''sudo apt-get install svn'')');
            else
                    url = sprintf('find\na version suitable for your system at\n\n\thttp://subversion.apache.org/packages.html');                    
            end
            fprintf(['\n Installation failed with error \n\n\t%s.\n\nIf this happened because you don''t have a working',...
                     ' copy of SVN, you can %s\n\n'],res,url)

        else
            addpath(grpath)
        end
        
    else
        stat = 1;
    end
    
    
else
    pth = fileparts(which('GazeReader'));
    fprintf('\nGazeReader installed at %s\n',pth)
    stat = 0;
end

gr_present = stat==0;
