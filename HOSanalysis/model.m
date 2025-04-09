
classdef model 
%
%  Regressor object for point-process models. 
% 
%  out=model(Response,RegIn,varargin)
%
%  Response is a struct with two fields:
%       .dat : Response data
%       .type: data type : 'binary' for binary data,'count' for count
%                          data,'times' for event times.
%       .fs : sampling rate.
% 
%  RegIn is a struct describing the independent variables.
%       .autodep: How to model the conditional dependence on the process history
%                 using a laguerre polynomial expansion. It should be a struct
%                 array with the following fields:
%                    .order: order of the expansion.
%                    .tau:   time constant.
%                 An empty array indicates no modeling of conditional dependence;
%                 with no dependence, statistical results are valid only if the
%                 process is Poisson. Default is
%                 .autodep = struct('order',{0 2},'tau',{0 .05 });

% C. Kovach 2019

    properties

       response;
       labels;
       modelType= 'binomial';
       autodep = struct('order',{0 4},'tau',{0 .025 });
       baseline= struct('order',4);
       event = struct('times',[],'evnt',[],'Trange',[]);
       sampling_rate=1;
       Tstart=0;
       timeBasis = 'polynomial';
       timeOrder = 16;
       regressors=regressor([]);
       fact2reg_args = {};
       intercept=true;
       codeincr = 0;
       do_llr_tests=false;
    end

    properties( Dependent = true)
       designMtx 

    end

    methods

        function me = model(Y,varargin)
          
            
            me.event(:) = [];
            
            if nargin < 1
                return
            elseif  isstruct(Y) && isfield(Y,'dat')
                me.response = Y.dat;
                me.sampling_rate = Y.fs;
                me.modelType = Y.type;
            elseif isa(Y,'model') || isstruct(Y) 
                fldn = setdiff(fieldnames(Y),{'designMtx'});
                for k = 1:length(fldn)
                    me.(fldn{k}) = Y.(fldn{k});
                end
            elseif isempty(Y)
                me = me([]);
                return
            else
                if isscalar(Y)
                    N = Y;
                    Y=false(N,1);
%                 else
%                     N = length(Y);
                end
            
                me.response = Y;          
                if isnumeric(varargin{1})
                    me.sampling_rate = varargin{1};
                    varargin(1)=[];
                end
            end
            i = 1;

            while i < length(varargin)

                switch(varargin{1})

                    case {'modelType','Tstart','autodep','event'}
                        out.(varargin{i}) = varargin{i+1};
                        i=i+1;
                    case {'eventTimes','eventTypes'}
                        fldn = lower(regexprep(varargin{i},'event',''));
                        me.event.(fldn) = varargin{i+1};
                    case 'Trange'
                         out.event.(varargin{i}) = varargin{i+1};
                    otherwise
                        if ~ischar(varargin{i})
                            error('Expecting a keyword string, found a %s object instead.',class(varargin{i}));
                        elseif ismember(varargin{i},fieldnames(me))
                            me.(varargin{i}) = varargin{i+1};
                            i = i+1;
                        else
                            error('%s is an unrecognized keyword.',varargin{i});
                        end
                end
                i=i+1;
            end
                
        end
        
        function out=get.designMtx(me)
            %%% Get the design matrix
            
            if ~isempty(me.regressors)
                codeincr = max([me.regressors.code]);
            else
                codeincr = me.codeincr;
            end

            %%% History regressor
            out = regressor([]);
            out(:) = [];
            Fs = zeros(size(me.response,1),length(me.event));
            for k = 1:length(me.event)
                
%                 evs = permute(me.event(k).evnt,[3 2 1]);
%                 evs = repmat(evs,size(evw.T,1),1);
                F = zeros(size(me.response,1),size(me.event(k).evnt,1));
                if ~isfield(me.event(k),'timeBasis')
                    me.event(k).timeBasis='polynomial';
                end
                switch me.event(k).timeBasis
                    case {'forward_laguerre','backward_laguerre'}
                        trts = round(me.event(k).times*me.sampling_rate)+1; %Zero corresponds to the 1st sample
                    otherwise
                        evw = me.get_event_window(me.event(k)); 
                        trts = evw.T(round(end/2),:);
                end
                for kk = 1:size(F,2)
                    [unq,~,unqi] = unique(me.event(k).evnt(kk,:));
                    F(trts,:) = unqi;
                end
%                 F(evw.T(:)*ones(1,size(evs,3))+size(F,1)*ones(numel(evw.T),1)*(0:size(evs,3)-1)) = evs(:);
                %Indicate times to be ignored outside the regression window
                mask = true(size(F,1),1);
%                 mask(evw.T) = false;
                  mask(trts,:) = false;
                if isfield(me.event(k),'center')
                    center = me.event(k).center;
                else
                    center = length(unique(me.event(k).evnt(kk,:)))>1;
                end
                if isfield(me.event(k),'label')
                    label = me.event(k).label;
                else
                    label = sprintf('Factor %i',k);
                end    
                if isfield(me.event(k),'fact2reg_args')
                    fact2reg_args = me.event(k).fact2regs_args;
                else
                    if ~isempty(out)
                        fact2reg_args={'codeincr',[out(end).code]};
                    else
                        fact2reg_args = {};
                    end
                end
                
                if ~iscell(label)
                    label = {label};
                end
%                 Fs(:,k) = F;
                
                
                Freg = fact2reg(F,'center',center,'ignore',mask,'labels',label,'fullintxn',fact2reg_args{:},'codeincr',codeincr);
                
                for kk = 1:length(Freg)
                    Freg(kk).value = fft(Freg(kk).value);
                end
                if center
                    trintcpt = zeros(size(Fs,1),1);
                    trintcpt(trts) = 1;
                    Freg(end+1) = regressor(fft(trintcpt),'label','trial intcpt','codeincr',Freg(end).code);
                end 
                Freg = Freg(arrayfun(@(x)~isempty(x.value),Freg));
                codeincr = Freg(end).code;
                %%%
                if ~isfield(me.event,'timeBasis')
                    me.event(k).timeBasis = me.timeBasis;
                end
                switch lower(me.event(k).timeBasis)
                    case {'forward_laguerre','backward_laguerre'}
                        TP = laguerreFilt((0:length(me.response)-1)'==0,me.event(k).laguerre_ord,me.event(k).laguerre_tau,me.sampling_rate);
                        TP = circshift(TP,[-1 0]);
                        Treg = regressor(fft(TP),'label','After','codeincr',codeincr);
                        if strcmpi(me.event(k).timeBasis,'backward_laguerre')
                          TP = flipud(TP);
                          Treg = regressor(fft(TP),'label','Before','codeincr',codeincr);
                        end

                    case {'polynomial','chebyt','bernstein'}
%                         timeOrder = [];
%                         if isfield(me.event(k),'timeOrder')
%                             timeOrder = me.event(k).timeOrder;
%                         end
%                         if isempty(timeOrder) %#ok<*PROP>
%                             timeOrder = me.timeOrder;
%                         end
%                         

                         TP = evw.P;
                         TP(size(me.response,1),:) = 0;
                         TP = circshift(TP,[-floor(size(evw.P,1)/2) 0]);


                              codeincr = max([codeincr,out.code]);
                        Treg = regressor(fft(TP),'label',sprintf('Window%i',k(k>1)),'codeincr',codeincr);
             
                    case 'boxcar'
                        error('Boxcar time basis is not implemented yet')
                    otherwise
                        error('%s is an unrecognized option for the time basis.')
                end
                FTreg=interaction(Freg,Treg);
                for kk = 1:length(Freg)
%                            FTreg(kk).window = k; 
                   FTreg(kk).window = me.event(k); 
                   FTreg(kk).value = ifft(FTreg(kk).value);
                end

                out = [FTreg,out]; %#ok<*AGROW>
            end   
            
            if ~isempty(me.autodep)  %%% Autoregressive component
                for k = 1:length(me.autodep)
%                     if isempty(out)
%                         codeincr=me.codeincr;
%                     else
                           codeincr = max([codeincr,out.code]);
%                     end

                     histreg = regressor(laguerreFilt(me.response,me.autodep(k),me.sampling_rate),'label',sprintf('historyLaguerre%i',k(k>1)),'codeincr',codeincr);
                    out(end+1) = histreg;
                end

            end
            if ~isempty(me.baseline)  %%% Autoregressive component
%                 if isempty(out)
%                     codeincr=me.codeincr;
%                 else
                    codeincr = max([codeincr,out.code]);
%                 end
                
                baseP = chebyT(length(me.response),me.baseline.order);
                basereg = regressor(baseP(:,2:end),'label',sprintf('Baseline'),'codeincr',codeincr);
                out(end+1) = basereg;

            end
   
           if ~isempty(me.regressors)
                out = [me.regressors,out];
           end
            for k = 1:length(out) %% Just force the codes to be unique and sequential
                out(k).code = k;
                out(k).codevec(:)=k;
            end

            me.codeincr=max([out.code]);
            %%% Factorial model
            
%             if me.intercept
%                 out(end+1) = regressor(ones(size(me.response,1),1),'label','intercept');
%             end
        end
        
        function evw = get_event_window(me,evnt)
            if nargin < 2 || isempty(evnt)
                evnt = me.event;
            elseif isnumeric(evnt)&& isscalar(evnt)
                evnt = me.event(evnt);
            end
            for k = 1:length(evnt)
                [evw(k).T,evw(k).tt] = chopper(evnt(k).Trange,evnt(k).times,me.sampling_rate);
                evw(k).T(evw(k).T<1)=1;
                evw(k).T(evw(k).T>length(me.response))=length(me.response);
                timeOrder = []; %#ok<*PROPLC>
                if isfield(evnt,'timeOrder')
                    timeOrder = evnt.timeOrder;
                end
                if isempty(timeOrder) %#ok<*PROP>
                    timeOrder = me.timeOrder;
                end
                if ~isfield(evnt,'timeBasis') || isempty(evnt(k).timeBasis)
                    evnt(k).timeBasis = me.timeBasis;
                end
                switch lower(evnt(k).timeBasis)
                    case 'bernstein'
                        evw(k).P = bernsteinp(length(evw(k).tt),timeOrder);
                    case {'chebyt','polynomial'}
                        evw(k).P = chebyT(length(evw(k).tt),timeOrder);
                    case {'forward_laguerre'}
                        evw(k).P = laguerreFilt(evw(k).tt(:)==0,evnt(k).laguerre_ord,evnt(k).laguerre_tau,me.sampling_rate);
                    case {'backward_laguerre'}
                        evw(k).P = flipud(laguerreFilt(flipud(evw(k).tt(:)==0),evnt(k).laguerre_ord,evnt(k).laguerre_tau,me.sampling_rate));
                end
%                 evw(k).intercept=evw(k).P(:,1);
%                 evw(k).P(:,1)=[];    
            end
        end
        
        function me=addregressor(me,varargin)
            
            dm = me.designMtx;
            me.regressors(end+1) = regressor(varargin{:},'codeincr',max([dm.code]));
            
        end
    end
end

