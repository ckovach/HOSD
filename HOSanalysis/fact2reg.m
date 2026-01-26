      
classdef fact2reg < regressor

% function RF  = fact2reg(F)
%
% Builds Regressor from a factor matrix, F.
%
% Each level in each factor generates a column in the regressor R.value
% which contains a dummy variable for the corresponding factor x level.
% RF is a regressor structure.
%
% Note that the columns of RF.value can be mapped onto the corresponding 
% factor x level interactions based on the colums of RF.factmat (factor) and
% RF.levmat (level). Which give the factors and levels for all terms in the
% interaction.                     
%
%
% function RF  = fact2reg(F,'option',value)
%
%  Create RF with the following options:
%         'intxnord'      -      Specify order of interaction (default  =1)
%         'fullintxn'     -      model with all interaction orders (takes
%                                 no argument)
%         'center'        -     Center regressors if true (default == true)
%                                   and decorrelated with the intercept.
%         'noptions'
%         'ignore'        - logical index into F: ignores entries in F for 
%                             which ignore is true. fact2reg does not create
%                             an indictor variable for unique entries
%                             contained only in ignored rows.
%
%   For remaining options see also REGRESSOR

% C Kovach 2008-2019
      
    properties 
       m
       ufw
       ufs
    end
    methods
        function me = fact2reg(F,varargin)
            
            me = me@regressor(0); 
            if isa(F,'regressor')
                me(1)=set(me(1),F(1));
                if length(F)>1
                    me = [me,fact2reg(F(2:end))];
                end
                return
            end
            intxnord = 1;
            center = true;
            noptions = [];
            labels = repmat({[]},1,size(F,2));
            % cellmeans = false;
            codeincr = 0;
            postmult = 1;
            i = 1;
            obsfreq = ones(size(F,1),1);

            ignore = false(size(F,1),1);

            varargin{end+1} = 'finis';
            while i <= length(varargin)
                switch lower(varargin{i})
                   case 'intxnord'                  %order of interaction
                        intxnord = varargin{i+1};
                        i = i+1;
                   case 'fullintxn'                 %mode with all interaction orders
                        intxnord = size(F,2);
                   case 'center'
                        center = varargin{i+1};
                        i = i+1;
                   case 'noptions'
                        noptions = varargin{i+1};
                        i = i+1;
                   case {'labels','label'}
                        labels = varargin{i+1};
                        i = i+1;
            %         case 'cellmeans'
            %             cellmeans = varargin{i+1};
            %             i = i+1;
                    case 'codeincr'
                        codeincr = varargin{i+1};
                        i = i+1;
                    case 'ignore'       % logical index into F: ignores entries in F for 
                                        % which ignore is true. fact2reg does not create
                                        % an indictor variable for unique entries
                                        % contained only in ignored rows.
                        ignore = varargin{i+1};
                        i = i+1;
                    case 'postmultiply'
                        postmult = varargin{i+1};
                        i = i+1;
                    case 'obsfreq'
                        obsfreq = varargin{i+1};
                        i = i+1;
                  case 'finis'

                  otherwise
                        error([varargin{i},' is not a valid keyword.'])
                end
                i = i+1;
            end 

            if ~iscell(labels)
                labels = {labels};
            end

            if iscell(F)  %Allow for cell arrays as well as numeric vectors
                eqfun = @(a,b) cellfun(@(a) isequal(a, b{1}),a);
            else
                eqfun = @(a,b) a==b;
            end

            for f = 1:size(F,2)

            %     me.ufs = unique(F(:,f));
                  me.ufs = unique(F(~ignore,f));  % Get unique values in the factor

                XF = zeros(size(F,1),length(me.ufs));

                for l = 1:length(me.ufs)

                    XF(:,l) = eqfun(F(:,f) , me.ufs(l) ).*~ignore;

                end


                if center
                    %remove mean
                    if ~all(obsfreq==1)
                %         u = ones(1,size(XF,1))'./sqrt(size(XF,1));
                        spbl = sparseblock(ones(1,size(XF,1)),noptions.*ones(size(obsfreq)));
                    else 
                        spbl = 1; 
                    end

                    u1 = (spbl'*obsfreq).*~ignore;
                    u1 = u1./sum(u1);
                    u2 = ~ignore;
                    m = u2*(u1'*XF);
                    XF = XF-m;
%                     XF = XF(:,2:end);
                    XF = XF(:,1:end-1);
                    me.m = m;
    %                 me.ufs = ufs;
    %                 fun= mkfun(ufs,m(find(~ignore,1),:));
    %             else
    %                 fun= mkfun(ufs);
    %                 me.ufs = ufs;
                end
                if isempty(XF)
                    continue
                end
                if isscalar(postmult)==1
                    XF = XF*postmult;
                else
                    XF = XF.*repmat(postmult,size(XF,1),1);
                end
                me(f) = set(me(1),regressor(XF,'noptions',noptions,'label',labels{f},'codeincr',codeincr));
%                 set(me(f),'value',XF,'noptions',noptions,'label',labels{f},'code',codeincr+1); %#ok<*AGROW>
                me(f).info.factorlabels = me.ufs;
%                 me(f).function =fun;
            end

            if intxnord>1
                me = cat(2,me,interaction(me,'intxnord',intxnord));
            end
        end
    %     %%%%%%%%%%%
    %     function fn = makefun(me)
    % 
    %    
    % 
    %         if isempty(me.m)
    %             fn = @(F)mkcfun(F,ufs);
    %         else %center
    %             fn = @(F)mkcfun(F,ufs,m);
    %         end
    % 
    %     end
        %%%%%%%%%
        function XF = makefun(me,F)

            if iscell(F)  %Allow for cell arrays as well as numeric vectors
                eqfun = @(a,b) cellfun(@(a) isequal(a, b{1}),a);
            else
                eqfun = @(a,b) a==b;
            end
            XF = zeros(size(F,1),length(me.ufs));

            for l = 1:length(me.ufs)

                XF(:,l) = eqfun(F , me.ufs(l) );

            end
            if ~isempty(me.m)
                XF = XF-repmat(me.m,size(XF,1),1);
                XF = XF(:,2:end);
            end

        end
        
       
    end

end

