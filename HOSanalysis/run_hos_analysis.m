
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
opts.save_space = false; %Save space by removing all HOS statistics, keeping only features and filters; 
opts.run_phase_randomized=false; %Run on phase randomized data as well as original for comparison
% opts.autodep = struct('order',8,'tau',.025);
outcode = char(java.util.UUID.randomUUID);
reseed;

% if false
    model;
% end

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
            outfn = sprintf('%s_%i_hos.mat',ld.blkdat.block,ld.chan(1).contact);
        otherwise  
            outfn = sprintf('%s_%i_out.mat',ld.blkdat.block,ld.chan(1).contact);
    end
    outputfile = fullfile(outputdir,outfn);
    
else
    useclust = false;
    inputdir='';
end
if nargin > 1 && isstruct(outputdir) || isfield(dat,'opts')
    if nargin > 1 && isstruct(outputdir)
        optsin = outputdir;
    else
        optsin = dat.opts;
    end
    
    fldn = fieldnames(optsin);
    for k = 1:length(fldn)
        opts.(fldn{k})=optsin.(fldn{k});
    end
    outputdir = '';
    outputfile = '';
end
if isempty(opts.resamp)
    [a,b] = rat(opts.target_fs/dat.fs(1),.01);
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
    thresh = iterz(dat.dat,6);
    [xdn,~,~,spk] = dbtDenoise(dat.dat,dat.fs(1),.1,'make plot',false,'spike window',.01);
    xdn =xdn +  0./(spk.filter>.5); % Place nans wherever the data are clipped
%     fr = getframe(fig);
%     dat.denoising = fr;
%     delete(fig);
%     shg
    dat.dat = xdn;
    dat.denoised=1;
end

dat.dat = double(dat.dat);
isn = isnan(dat.dat);
dat.dat(isn) = 0;
% fprintf('\ndat.dat type is: %s, size:[%i %i]',class(dat.dat),size(dat.dat));
dat.dat = resample(dat.dat,opts.resamp(1),opts.resamp(2)) + 0./(resample(double(isn),opts.resamp(1),opts.resamp(2))==0);
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
%    z = zscore(dat.dat);
%    discard = any(z(T)>opts.zthresh);
    z = iterz(dat.dat,opts.zthresh);
   discard = any(isnan(z(T)));
   if mean(discard)<.5 % Run on all data if too many windows are discarded.
       segment.wint(discard)=[];
   end
end

hos(1) = hosobject(opts.hos_order);
hos(opts.ncomp) = hosobject(opts.hos_order);

fldn = fieldnames(opts);
fldn = fldn(ismember(fldn,fieldnames(hos)));
for k = 1:length(fldn)
    opts.hosargs(end+1:end+2) = [fldn(k),{opts.(fldn{k})}];
end
seed = reseed;

z = zscore(double(dat.dat));
hosargs = [{size(T,1),dat.fs(1),opts.lowpass,[],[]},opts.hosargs];
if nargin > 1 && exist('outputfile','var')&&exist(outputfile,'file') && ~opts.redo_hosd && ~opts.save_space
    try
        load(outputfile,'hos','segment')
    catch err
        delete(outputfile)
        error(err)
    end
    segment = segment(1);
    segment.wintadj=[];
else
     hos.initialize(hosargs{:});
%     if apply_to_chopped_data
%         hos.get_block(z(T));  % Deflation is done on the chopped data
%       
%     else
        hos.get_block(z,[],[],segment);  % This allows the entire record to be used in the deflation step.
      
%     end
    if opts.run_phase_randomized
        hosphaserand = hosobject(hos);
        znz = z;
        znz(isnan(z))=0;
        zrand = zscore(real(ifft(exp(2*pi*1i*rand(size(z))).*abs(fft(znz)))));
        zrand(isnan(z)) = nan;
        
        if apply_to_chopped_data
            hosphaserand.get_block(zrand(T)); 
        else
            hosphaserand.get_block(zrand,[],[],segment);
        end
    end
end
% xrec = hos.xrec(z); % Get the reconstruction
% ximp = hos.ximp(z);
% xfilt = hos.xfilt(z);

bsidout(1).hos= hos;


if opts.run_phase_randomized
    bsidout(1).hosphaserand= hosminimal(hosphaserand);
end
bsidout(1).hosargs = [{opts.hos_order},hosargs];
bsidout(1).dat = z;
if ~isfield(dat,'chan')
    dat.chan = struct('label','Channel','number',dat.ChannelNumber(1)+1,'channel',dat.ChannelNumber(1)+1,'contact',dat.ChannelNumber(1)+1);
end
if ~isfield(dat,'block')
    dat.block = struct('block','???-???','subprotocol','????');
end

bsidout(1).block = dat.block;
bsidout(1).chan = dat.chan;
bsidout(1).fs = dat.fs(1);
% bsidout(1).origdatafile=dat.origdatafile;
bsidout(1).opts = opts;
bsidout(1).segment = segment;
bsidout(1).randseed = seed;
if exist(inputdir,'dir')&& (~isfield(opts,'no_anls') || ~opts.no_anls)
%     res = hos_regression_analysis(bsidout);
    res = feval(opts.stats_function,bsidout,inputdir);
else
    res=[];
end

if isfield(opts,'hos_regressor')&& ~isempty(opts.hos_regressor)
        res.hosregresult = hos.hos_regress(bsidin.dat,opts.hos_regressor);
end

bsidout.result = res;
bsidout(1).dat = z;
    

    
% bsidout.hos= hos;
% bsidout(1).dat = z;
% bsidout(1).chan = dat.chan;
% bsidout(1).block = dat.block;
% bsidout(1).fs = dat.fs(1);
% % bsidout(1).origdatafile=dat.origdatafile;
% bsidout(1).opts = opts;


if isfield(opts,'make_plots') && opts.make_plots && nargin >1
    if ~isfield(opts,'plot_function')
       bsidout(1).res = hos_regression_plot(bsidout,useclust,outputdir);
    else
        bsidout(1).res = feval(opts.plot_function,bsidout,useclust,outputdir);
    end
end
if opts.save_space
    bsidout(1).hos= hosminimal(hos);
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
    
