function res = hos_regression_analysis(bsidin,varargin)

opts = bsidin.opts;
hos = bsidin.hos;
segment = bsidin.segment(1);
xthresh = bsidin.hos.xthresh(bsidin.dat).^2;

res = struct('atts',[],'Mbi',[],'afrqs',[],'segment',[],'fit',[],'model',[]);

if isfield(opts,'modelopts') && ~isempty(opts.modelopts)
    mdl = model(opts.modelopts);
   mdl.event = opts.modelopts.event;
    mdl.sampling_rate=bsidin.fs(1);

elseif isfield(opts,'event')
        
    mdl = model;
    mdl.event = opts.event;
    mdl.sampling_rate=bsidin.fs(1);
else
    mdl = [];
end

if ~isfield(opts,'no_anls') || ~opts.no_anls
    x = bsidin.dat;
    x(isnan(x))=0;
    for compi = 1:length(hos)
        segment = hos(compi).segment;
        if isempty(hos(compi).segment.wint)
            if isfield(opts,'windur') && isstruct(opts.windur)
                hos(compi).segment = opts.windur;
            end
            if isempty(hos(compi).segment.wint)        
                continue
            end
        end
%         segment.wintadj = hos(compi).delay + segment.wint;
        for bi = 1:size(opts.bands,1)   

           dbx = dbt(x,bsidin.fs(1),opts.bands(bi,3),'upsample',4,'lowpass',opts.bands(bi,2),'highpass',opts.bands(bi,1),'remodphase',true,'centerDC',false);
            %%% envelope smoothing

            dbx.blrep = dbx.blrep./abs(dbx.blrep).*sqrt(convn(abs(dbx.blrep).^2,hann(3*opts.time_freq_smoothn),'same'));
           [As{bi},attsb] = choptf(segment.Trange*segment.fs/bsidin.fs(1),segment.wintadj*segment.fs/bsidin.fs(1),dbx,segment.Trange*segment.fs/bsidin.fs(1)); %#ok<*AGROW>
           atts{bi} = attsb-mean(segment.Trange(:)*segment.fs/bsidin.fs(1));
           Mbi{bi} = mean(20*log10(abs(As{bi})),3)';
        %    Mevbi{bi} = 20*log10(abs(mean(As{bi},3)))';
           frqs{bi}=dbx.frequency;
        end           
        [~,pks] = getpeak2(xthresh(:,compi));
         imp = full(pks==1);

    %      opts.autodep = struct('order',{0 8},'tau',{0 , median(diff(find(imp)))/bsidin.fs});

        if isfield(opts,'inpt')
           inpt = opts.inpt;
           dbinpt = dbt(inpt.dat,inpt.fs,20,'lowpass',min(inpt.fs/2,4e3));
    %        imp = full(ximp(:,compi));

           dbsnd = dbt(abs(dbinpt.blrep),dbinpt.sampling_rate,.25);
%            sndY = reshape(dbsnd.blrep,length(dbsnd.time),numel(dbsnd.blrep(1,:,:)));
           dbimp = dbt(imp,bsidin.fs(1),.25,'lowpass', dbsnd.lowpass);
%            impX = reshape(dbimp.blrep,length(dbimp.time),numel(dbimp.blrep(1,:,:)));

           coh = dbtcoh(dbimp,dbsnd);
           C = squeeze(coh);
           C(:,2*end+1) = 0;
           CTF = fftshift(real(ifft(C,[],2)),2);
           res(compi).CTF = CTF;

           res(compi).ctft = ((0:size(C,2)-1)-floor(size(C,2)/2))./size(C,2)./diff(dbimp.frequency(1:2));
            res(compi).ctftfrq = dbinpt.frequency';
         end


        res(compi).atts=atts;
        res(compi).Mbi=Mbi;
        res(compi).afrqs = frqs;

        if ~isempty(mdl)
             opts.autodep = struct('order',{ 8},'tau',{ median(diff(find(imp)))/bsidin.fs(1)});
             mdl.autodep = opts.autodep;

            mdl.response = imp;

            if isfield(opts,'regressors')
                mdl.addregressor(opts.regressors);
            end

            if ~all(isnan(imp)) &&( ~isempty(mdl.regressors) || any(imp(mdl.get_event_window(1).T(:))))
                fit = fitmod(mdl);            
                res(compi).fit = fit;
            else
                res(compi).fit = [];
            end            
                res(compi).model = model;
        end
       res(compi).segment = segment;
    end
%     if isfield(opts,'hos_regressor')&& ~isempty(opts.hos_regressor)
%             res.hosregresult = hos.hos_regress(bsidin.dat,opts.hos_regressor);
%     end
else
    res=[];
end