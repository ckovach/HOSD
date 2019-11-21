




function bsidout=run_hos_analysis(dat, outputdir,inputfiles,jobindex)

opts.lowpass = 200;
opts.windur = 3;
opts.povlp=.5;
opts.target_fs = 500;
opts.resamp = [];
opts.bsidargs={};
opts.bands = [150 200 20
        80 125  15
        40 80   10
        20 40   5
        10 20   2
        0  10   1];
opts.time_freq_smoothn=5;
opts.ncomp = 2;    
opts.zthresh = Inf; %Exclude windows in which max value exceeds this z-score threshold.
% opts.do_regression = false;
opts.version = 'hos';
opts.hos_order=3;
opts.nperm = 5e3;
opts.hosargs = {};
opts.make_plots = true;
opts.redo_hosd = false;
opts.dbt_denoise = true;
% opts.autodep = struct('order',8,'tau',.025);
outcode = char(java.util.UUID.randomUUID);
reseed;

sl = license('checkout','signal_toolbox');
atn=1;
while isequal(sl,0)
    fprintf('\nFailed to checout the signal processing toolbox. Attempt %i',atn);
    atn=atn+1;
    
    sl = license('checkout','signal_toolbox');
    pause(1)
end

t0=tic;
if ischar(dat)
    inputdir = dat;

    if jobindex==0; jobindex =1 ;end
    useclust=true;
    
    [~,fn,ext] = fileparts(inputfiles{jobindex});
    switch ext
        case '.mat'
            ld = load(fullfile(inputdir,fn));
            if isnumeric(ld.dat)
                dat = ld;
            else
                dat = ld.dat;
            end
            if isfield(ld,'opts')
                fldn = fieldnames(ld.opts);
                for k = 1:length(fldn)
                    opts.(fldn{k})=ld.opts.(fldn{k});
                end
            end
            if isfield(ld,'chan')
                dat.chan = ld.chan;
            end
            if isfield(ld,'blkdat')
                dat.block = ld.blkdat;
            end
             if isfield(ld,'block')
                dat.block = ld.block;
             end
            if ~isfield(dat.block,'subprotocol')
                dat.block.subprotocol='';
            end
        case '.ncs'
            dat = readncs([fn,ext],inputdir);
    end
    fid = fopen(fullfile(outputdir,'manifest.txt'),'a+');
    fprintf(fid,'\n%s\t0\tINPUT\t%s\t%0.3fs',[fn,ext],outcode,toc(t0));
    fclose(fid);
    if exist(fullfile(inputdir,'model.mat'),'file')
        ldopt=load(fullfile(inputdir,'model.mat'));
        fldn = fieldnames(ldopt.opts);
        for k = 1:length(fldn)
            opts.(fldn{k})=ldopt.opts.(fldn{k});
        end
        if strcmp(ext,'.ncs')
             chn = regexp(dat.file,'LFPx(\d*)','tokens','once');
            chn = str2double(chn{1});
            chan = opts.block.lozchannels([opts.block.lozchannels.channel]==chn);
            ld.chan = chan;
            ld.blkdat = opts.block;
            dat.chan = chan;
            dat.block = opts.block;
        elseif ~isfield(ld,'blkdat') && isfield(ld,'block')
            ld.blkdat=ld.block;
        end
            
        opts.modelopts = ldopt.model;
    end
    
    
    switch opts.version
        case 'hos'
            outfn = sprintf('%s_%i_hos.mat',ld.blkdat.block,ld.chan.contact);
        otherwise  
            outfn = sprintf('%s_%i_out.mat',ld.blkdat.block,ld.chan.contact);
    end
    outputfile = fullfile(outputdir,outfn);
    
else
    useclust = false;
end
if nargin > 1 && isstruct(outputdir)
    optsin = outputdir;
    fldn = fieldnames(optsin);
    for k = 1:length(fldn)
        opts.(fldn{k})=optsin.(fldn{k});
    end
    outputdir = '';
    outputfile = '';
end
if isempty(opts.resamp)
    [a,b] = rat(opts.target_fs/dat.fs(1),.1);
    if b>a
        opts.resamp = [a b];
    else
        opts.resamp = [1 1];
    end
elseif length(opts.resamp)==1  % Scalar value for resampling is treated as decimation factor
    opts.resamp = [1 opts.resamp];
end

if opts.dbt_denoise && (~isfield(dat,'denoised')  ||  ~dat.denoised)
%     fig = figure;
    xdn = dbtDenoise(dat.dat,dat.fs(1),.1,'make plot',false,'spike window',.01);
%     fr = getframe(fig);
%     dat.denoising = fr;
%     delete(fig);
%     shg
    dat.dat = xdn;
    dat.denoised=1;
end

dat.dat = resample(double(dat.dat),opts.resamp(1),opts.resamp(2));
dat.fs = dat.fs(1)*opts.resamp(1)./opts.resamp(2);

n = length(dat.dat);
% 
% if isfield(opts,'modelopts') && ~isempty(opts.modelopts)
%     mdl = model(opts.modelopts);
%    mdl.event = opts.modelopts.event;
%     mdl.sampling_rate=dat.fs(1);
% 
% elseif isfield(opts,'event')
%         
%     mdl = model;
%     mdl.event = opts.event;
%     mdl.sampling_rate=dat.fs(1);
% else
%     mdl = [];
% end


apply_to_chopped_data = isstruct(opts.windur);

if ~apply_to_chopped_data
    segment = dat.fs(1)*opts.windur;
    
    Trange= [-1 1]*segment/2;
    
    segment = struct('Trange',Trange,'fs',1,'povlp',opts.povlp);
    segment.wint= 1/segment.fs:diff(segment.Trange)*(1-segment.povlp):n/segment.fs;
    
else
    segment = opts.windur;
end

[T,~] = chopper(segment.Trange,segment.wint,segment.fs);
T(T<1)=1;T(T>n)=n;

zscore = @(x)(x-nanmean(x))./nanstd(x);

segment.wint(any(isnan(dat.dat(T))))=[];
T(:,any(isnan(dat.dat(T))))=[];

if opts.zthresh<Inf
   z = zscore(dat.dat);
   discard = any(z(T)>opts.zthresh);
   segment.wint(discard)=[];
end


hos(1) = hosobject(opts.hos_order);
hos(opts.ncomp) = hosobject(opts.hos_order);
z = zscore(double(dat.dat));
if nargin > 1 && exist('outputfile','var')&&exist(outputfile,'file') && ~opts.redo_hosd
    load(outputfile,'hos','segment')
    segment = segment(1);
    segment.wintadj=[];
else
     hos.initialize(size(T,1),dat.fs(1),opts.lowpass,[],[],opts.hosargs{:});
    if apply_to_chopped_data
        hos.get_block(z(T));  % Deflation is done on the chopped data
    else
        hos.get_block(z,[],[],segment);  % This allows the entire record to be used in the deflation step.
    end
end
% xrec = hos.xrec(z); % Get the reconstruction
% ximp = hos.ximp(z);
% xfilt = hos.xfilt(z);

bsidout.hos= hos;
bsidout(1).dat = z;
bsidout(1).chan = dat.chan;
bsidout(1).block = dat.block;
bsidout(1).fs = dat.fs(1);
% bsidout(1).origdatafile=dat.origdatafile;
bsidout(1).opts = opts;
bsidout(1).segment = segment;

if ~isfield(opts,'no_anls') || ~opts.no_anls
%     res = hos_regression_analysis(bsidout);
    res = feval(opts.stats_function,bsidout);
end

if isfield(opts,'hos_regressor')&& ~isempty(opts.hos_regressor)
        res.hosregresult = hos.hos_regress(bsidin.dat,opts.hos_regressor);
end
% xthresh = hos.xthresh(z);
% 
% if ~isfield(opts,'no_anls') || ~opts.no_anls
%     x = dat.dat;
%     x(isnan(x))=0;
%     for compi = 1:length(hos)
%         segment.wintadj = hos(compi).delay + segment.wint;
%         for bi = 1:size(opts.bands,1)   
% 
%            dbx = dbt(x,dat.fs(1),opts.bands(bi,3),'upsample',4,'lowpass',opts.bands(bi,2),'highpass',opts.bands(bi,1),'remodphase',true,'centerDC',false);
%             %%% envelope smoothing
% 
%             dbx.blrep = dbx.blrep./abs(dbx.blrep).*sqrt(convn(abs(dbx.blrep).^2,hann(3*opts.time_freq_smoothn),'same'));
%            [As{bi},attsb] = choptf(segment.Trange*segment.fs/dat.fs(1),segment.wintadj*segment.fs/dat.fs(1),dbx,segment.Trange*segment.fs/dat.fs(1)); %#ok<*AGROW>
%            atts{bi} = attsb-mean(segment.Trange(:)*segment.fs/dat.fs(1));
%            Mbi{bi} = mean(20*log10(abs(As{bi})),3)';
%         %    Mevbi{bi} = 20*log10(abs(mean(As{bi},3)))';
%            frqs{bi}=dbx.frequency;
%         end           
%         [~,pks] = getpeak2(xthresh(:,compi));
%          imp = full(pks==1);
% 
%     %      opts.autodep = struct('order',{0 8},'tau',{0 , median(diff(find(imp)))/dat.fs});
% 
%         if isfield(opts,'inpt')
%            inpt = opts.inpt;
%            dbinpt = dbt(inpt.dat,inpt.fs,20,'lowpass',min(inpt.fs/2,4e3));
%     %        imp = full(ximp(:,compi));
% 
%            dbsnd = dbt(abs(dbinpt.blrep),dbinpt.sampling_rate,.25);
% %            sndY = reshape(dbsnd.blrep,length(dbsnd.time),numel(dbsnd.blrep(1,:,:)));
%            dbimp = dbt(imp,dat.fs(1),.25,'lowpass', dbsnd.lowpass);
% %            impX = reshape(dbimp.blrep,length(dbimp.time),numel(dbimp.blrep(1,:,:)));
% 
%            coh = dbtcoh(dbimp,dbsnd);
%            C = squeeze(coh);
%            C(:,2*end+1) = 0;
%            CTF = fftshift(real(ifft(C,[],2)),2);
%            res.CTF = CTF;
% 
%            res.ctft = ((0:size(C,2)-1)-floor(size(C,2)/2))./size(C,2)./diff(dbimp.frequency(1:2));
%             res.ctftfrq = dbinpt.frequency';
%          end
% 
% 
%         res(compi).atts=atts;
%         res(compi).Mbi=Mbi;
%         res(compi).afrqs = frqs;
% 
%         if ~isempty(mdl)
%              opts.autodep = struct('order',{ 8},'tau',{ median(diff(find(imp)))/dat.fs(1)});
%              mdl.autodep = opts.autodep;
% 
%             mdl.response = imp;
% 
%             if isfield(opts,'regressors')
%                 mdl.addregressor(opts.regressors);
%             end
% 
%             if ~all(isnan(imp)) &&( ~isempty(mdl.regressors) || any(imp(mdl.get_event_window.T(:))))
%                 fit = fitmod(mdl);            
%                 res(compi).fit = fit;
%             else
%                 res(compi).fit = [];
%             end            
%                 res(compi).model = model;
%         end
%         bsidout(1).segment(compi) = segment;
%     end
%     if isfield(opts,'hos_regressor')&& ~isempty(opts.hos_regressor)
%             bsidout.hosregresult = hos.hos_regress(z,opts.hos_regressor);
%     end
% else
%     res=[];
% end
% bsidout.hos= hos;
bsidout.result = res;
bsidout(1).dat = z;
    

    
% bsidout.hos= hos;
% bsidout(1).dat = z;
% bsidout(1).chan = dat.chan;
% bsidout(1).block = dat.block;
% bsidout(1).fs = dat.fs(1);
% % bsidout(1).origdatafile=dat.origdatafile;
% bsidout(1).opts = opts;


if isfield(opts,'make_plots') && opts.make_plots
    if ~isfield(opts,'plot_function')
       bsidout(1).res = hos_regression_plot(bsidout,useclust,outputdir);
    else
        bsidout(1).res = feval(opts.plot_function,bsidout,useclust,outputdir);
    end
end

try
    fid = fopen([mfilename,'.m']);
    bsidout(1).COM = fread(fid,'uchar=>char')';
    fclose(fid);
catch
    bsidout(1).COM = '';
end
if useclust
    bsidout = stripfunctions(bsidout);
    save(outputfile,'-struct','bsidout');
    fid = fopen(fullfile(outputdir,'manifest.txt'),'a+');
    fprintf(fid,'\n%s\t0\tOUTPUT\t%s\t%0.3fs',outfn,outcode,toc(t0));
    fclose(fid);
    

end
    
