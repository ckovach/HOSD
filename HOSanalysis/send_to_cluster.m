function xne = send_to_cluster(files, model,optsin)

% xne = send_to_cluster(files, model,opts)
% 
% Send a batch of run_hos_analysis jobs to the cluster. 
%
% files - input data as separate files with 1 channel per file. 
%         These can be neuralynx files.
% model - is a model object containing a prototype model
%
% opts  - structure with xargon options. If opts is a string, it will be 
%         assumed to be the output directory.
%
% See RUN_HOS_ANALYSIS, MODEL, REGRESSOR

opts.queue = 'UI,CCOM,all.q';
opts.profile = 'mid_mem';
opts.skipdone = true;
opts.nslots = 4;
opts.concatenate = [];
opts.window = 'hann';

if nargin < 2
    model = [];
end
if isa(files,'xargon') || isa(files,'xneon')
    xne = files;
    files = {xne.datafiles(1:end-1).orig};
%     mdlfile = xne.datafiles(end).orig;
%     mdlfile = fullfile(xne.datafiles(end).orig;
    if ~isempty(xne.subpaths.assets.local)
        mdlfile = fullfile(xne.subpaths.assets.local,'model.mat');
    end
    if ~exist(mdlfile,'file')        
        mdlfile = fullfile(xne.local_save_dir,'options.mat');
    end
    if ~exist(mdlfile,'file')        
        mdlfile = fullfile(xne.local_save_dir,'model.mat');
    end
    ld = load(mdlfile,'opts','model');
    if nargin < 3 || isempty(optsin)
        optsin = ld.opts;
    end
    if nargin < 2 || isempty(model)
        model = ld.model;
    end
else
    if ischar(files)
        [~,fn,ext] = fileparts(files);
        if isempty(ext)
            files = get_protocol_labwiki(files);
        else
            files = {files};
        end
    end

    if isstruct(files)
       block = locateNlx(files);
       files = fullfile(block.blkfiles.path,block.blkfiles.lfp);
       opts.block = block;
    end

    xne = xargon(which('run_hos_analysis'));
end
if ~isempty(model) 
    if nargin > 2 && isstruct(optsin)
        fldn = fieldnames(optsin);
        for k = 1:length(fldn)
            opts.(fldn{k}) = optsin.(fldn{k});
        end
    elseif ischar(optsin)
        opts.savedir=optsin;
        xne.local_save_dir = opts.savedir;
    end
end

fldn = fieldnames(opts);
for k = 1:length(fldn)
    if ismember(fldn{k},fieldnames(xne))
      xne.(fldn{k}) = opts.(fldn{k});
      opts = rmfield(opts,fldn{k});
    end
end
if opts.skipdone && exist(xne.subpaths.output.local,'dir') 
    existing_dir = xne.subpaths.output.local;
elseif opts.skipdone && exist(xne.local_save_dir,'dir') 
    existing_dir = xne.local_save_dir;
else
    existing_dir = '';
end


if ~isempty(existing_dir) && exist(fullfile(existing_dir,'manifest.txt'),'file')
    manfile = fullfile(existing_dir,'manifest.txt');
    fid = fopen(manfile,'r');
    txt = fread(fid,'uchar=>char')';
    fclose(fid);
    infiles = regexp(txt,'\n([\w_.\-]*)\t*0\t*INPUT\t*([\w-]*)','tokens');
    infiles = cat(1,infiles{:});
    outfiles = regexp(txt,'\n([\w_.\-]*)\t*0\t*OUTPUT\t*([\w-]*)','tokens');
    outfiles = cat(1,outfiles{:});
    if isempty(outfiles)
                missing = true(size(files));
    else

        missing_files = cellfun(@(x)~exist(fullfile(existing_dir,x),'file')&~exist(fullfile(existing_dir,'figs',x),'file'),outfiles(:,1));
        if exist('xne','var') && isa(xne,'xargon') && exist(xne.local_save_dir,'dir')
            missing_files = missing_files & cellfun(@(x)~exist(fullfile(xne.local_save_dir,x),'file')&~exist(fullfile(xne.local_save_dir,'figs',x),'file'),outfiles(:,1));        
        end

        outfiles(missing_files,:)={'xxxx'};
    %     donefiles = infiles(ismember(infiles(:,2),outfiles(:,2)),1);
        if ~isempty(infiles)||isempty(outfiles)
        [ism,ismi]= ismember(infiles(:,2),outfiles(:,2));
        if any(ism)
            ddat = dir(fullfile(existing_dir,'*hos.mat'));
            lddat = load(fullfile(existing_dir,ddat(1).name));
            donefiles = infiles(ism);
            pdffiles = outfiles(contains(outfiles(:,1),'pdf'));
            pdfco = regexp(pdffiles,'contact_(\d*)_bispectral','tokens','once');
            pdfco = cellfun(@str2num,[pdfco{:}]);
            [~,c2ci] = ismember(pdfco,[lddat.block.lozchannels.contact]);
            pdfch = [lddat.block.lozchannels(c2ci(c2ci~=0)).channel];
            [~,ff,ext] = cellfun(@fileparts,files,'uniformoutput',false);
            [fism,ford] = ismember(strcat(ff,ext),infiles(ism,1));
            donech = regexp(outfiles(:,1),'_(\d*)_hos','tokens','once');
             [donech,~,unqi] = unique(cellfun(@str2double,[donech{ismi(ism)}]),'stable');
    %         donech = cellfun(@str2double,[donech{ismi(ism)}]);
        %     files = files(~ismember(strcat(ff,ext),donefiles));
            missing = ~ismember(strcat(ff,ext),donefiles);
        else
            missing = true(size(files));
        end
        else
            missing = true(size(files));
        end        
        if any(missing)
            xne.jobindices=find( missing );
        else
           donech(~missing) =donech(unqi(ford(~missing)));

            xne.jobindices = find(~isempty(pdfch) & ~ismember(donech,pdfch));
        end

    %    donef = dir(fullfile(existing_dir,'*_hos.mat'));
    %    donech = regexp({donef.name},'_(\d*)_hos[.]mat','tokens','once');
    %    donech = cellfun(@str2double,donech);
    % %    chs = [opts.block.lozchannels,opts.block.hizchannels];
    %    [~,fns] = cellfun(@fileparts,files,'uniformoutput',false);
    %    availch = regexp(fns,'(\d*)$','tokens','once');
    %    availch = cellfun(@str2double,[availch{:}]);
    %    undoneidx = find(~ismember(availch,donech));
       if isempty(xne.jobindices)
           fprintf('\nAll channels finished')
           if ~strcmp(xne.jobindices,'none')
               xne.finish();
           end
           return
       elseif xne.nparallel <length(infiles)
           fprintf('\n%i channels remaining of %i in %s',xne.nparallel,length(infiles),existing_dir);
       end
    end
%     xne.jobindices=undoneidx;  
else
    manfile = '';
end
if isempty(xne.jobindices)
    xne.nparallel = length(files);
end
% 
% fldn = fieldnames(opts);
% for k = 1:length(fldn)
%     if ismember(fldn{k},fieldnames(xne))
%       xne.(fldn{k}) = opts.(fldn{k});
%       opts = rmfield(opts,fldn{k});
%     end
% end

if nargin > 1 || ~isempty(model)
    mdlfile = fullfile(xne.tempdir,'model.mat');
    save(mdlfile,'model','opts');
    xne.datafiles = [files,{mdlfile}];
else
    xne.datafiles = files;

end

xne.rerun_if_aborted='yes';

xne.create_job;
xne.make_bash_script;
xne.make_matlab_wrapper;

if ~isempty(manfile)
    try
        copyfile(manfile,xne.subpaths.output.local);
    end
end

xne.finish = @(varargin)finish(xne,varargin{:});

function finish(xne,varargin)

optsfile = fullfile(xne.subpaths.assets.local,'model.mat');
if exist(optsfile)
    if ~exist(xne.local_save_dir,'dir')
        mkdir(xne.local_save_dir)
    end
    copyfile(optsfile,xne.local_save_dir);
end

xne.default_finish();

manfile = fullfile(xne.local_save_dir,'manifest.txt');
if exist(manfile,'file')
    fid = fopen(manfile,'r');
    txt = fread(fid,'uchar=>char')';
    fclose(fid);
    re = regexp(txt,'([^\n\s]*[.]pdf)[^\n]','tokens');
    re =[re{:}];
    re = unique(re);
    re2 = regexp(re,'contact_(\d*)_','tokens','once');
    cnum = cellfun(@(x)str2num(['0',x{:}]),re2);
    [srt,srti] = sort(cnum);
        
    fns = fullfile(xne.local_save_dir,'figs',re(srti(srt>0)));
    com = sprintf('gs -sDEVICE=pdfwrite -dEPSCrop  -dMaxInlineImageSize=100000 -o%s%s%s_summary.pdf %s',xne.local_save_dir,filesep,regexp(xne.local_save_dir,['[^',filesep,']*$'],'match','once'),sprintf(' %s ',fns{:}));
    [err,out]=system(com);
 
    fprintf('%s',out)
    if err~=0
        fprintf('\n\nThe ghostscript command to compile figures into a pdf appears to have failed.\nSee above for clues.\nThe command was as follows:\n\n\t%s',com);
    end
     
end

     
     