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
opts.nslots = 4;
opts.skipdone = true;

if isa(files,'xargon')
    xne = files;
    files = {xne.datafiles(1:end-1).orig};
    mdlfile = xne.datafiles(end).orig;

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
if nargin >2 
    if isstruct(optsin)
        fldn = fieldnames(optsin);
        for k = 1:length(fldn)
            opts.(fldn{k}) = optsin.(fldn{k});
        end
    else
        opts.savedir=optsin;
        xne.local_save_dir = opts.savedir;
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
    infiles = regexp(txt,'\n([\w_-.]*)\t0\tINPUT\t([\w-]*)','tokens');
    infiles = cat(1,infiles{:});
    outfiles = regexp(txt,'\n([\w_-.]*)\t0\tOUTPUT\t([\w-]*)','tokens');
    outfiles = cat(1,outfiles{:});
    donefiles = infiles(ismember(infiles(:,2),outfiles(:,2)),1);
    files = files(~ismember(files,donefiles));
    
%    donef = dir(fullfile(existing_dir,'*_hos.mat'));
%    donech = regexp({donef.name},'_(\d*)_hos[.]mat','tokens','once');
%    donech = cellfun(@str2double,donech);
% %    chs = [opts.block.lozchannels,opts.block.hizchannels];
%    [~,fns] = cellfun(@fileparts,files,'uniformoutput',false);
%    availch = regexp(fns,'(\d*)$','tokens','once');
%    availch = cellfun(@str2double,[availch{:}]);
%    undoneidx = find(~ismember(availch,donech));
   if isempty(files)
       fprintf('\nAll channels finished')
       return
   elseif length(files) <length(infiles)
       fprintf('\n%i channels remaining of %i',length(files),length(infiles));
   end
%     xne.jobindices=undoneidx;  
else
    manfile = '';
end

xne.nparallel = length(files);

fldn = fieldnames(opts);
for k = 1:length(fldn)
    if ismember(fldn{k},fieldnames(xne))
      xne.(fldn{k}) = opts.(fldn{k});
      opts = rmfield(opts,fldn{k});
    end
end

if nargin > 1 
    mdlfile = fullfile(xne.tempdir,'model.mat');
    save(mdlfile,'model','opts');
end

xne.rerun_if_aborted='yes';

xne.datafiles = [files,{mdlfile}];
xne.create_job;

if ~isempty(manfile)
    copyfile(manfile,xne.subpaths.output.local);
end

xne.finish = @(varargin)finish(xne,varargin{:});

function finish(xne,varargin)

xne.default_finish();

manfile = fullfile(xne.local_save_dir,'manifest.txt');
if exist(manfile,'file')
    fid = fopen(manfile,'r');
    txt = fread(fid,'uchar=>char')';
    re = regexp(txt,'([^\n\s]*[.]pdf)[^\n]','tokens');
    re =[re{:}];
    re2 = regexp(re,'contact_(\d*)_','tokens','once');
    cnum = cellfun(@(x)str2num(['0',x{1}]),re2);
    [srt,srti] = sort(cnum);
        
    fns = fullfile(xne.local_save_dir,'figs',re(srti));
    com = sprintf('gs -sDEVICE=pdfwrite -dEPSCrop  -dMaxInlineImageSize=100000 -o%s%s%s_summary.pdf %s',xne.local_save_dir,filesep,xne.local_save_dir,sprintf(' %s ',fns{:}));
    [err,out]=system(com);
 
    fprintf('%s',out)
    if err~=0
        fprintf('\n\nThe ghostscript command to compile figures into a pdf appears to have failed.\nSee above for clues.\nThe command was as follows:\n\n\t%s',com);
    end
     
end

     
     