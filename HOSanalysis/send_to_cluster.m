
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
opts.plot_function = 'hos_regression_plot';
opts.stats_function = 'hos_regression_analysis';

if nargin < 2
    model = [];
end
if  isa(files,'xargon') || isa(files,'xneon')
    xne = files;
    opts.queue = xne.queue;
    opts.profile = xne.profile;
    opts.nslots = xne.nslots;
    
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
    if exist(mdlfile,'file')
        ld = load(mdlfile,'opts','model');
  
        if nargin < 3 || isempty(optsin)
            optsin = ld.opts;
        end
        if nargin < 2 || isempty(model)
            model = ld.model;
        end
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
       files = {};
       for bk = 1:length(block.blkfiles)
           files = [files,fullfile(block.blkfiles(bk).path,block.blkfiles(bk).lfp)];
       end
       opts.block = block;
    end

    xne = xargon(which('run_hos_analysis'));
end
if exist('optsin','var') &&~isempty(optsin) 
    if  isstruct(optsin)
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


if ~isempty(existing_dir) % && exist(fullfile(existing_dir,'manifest.txt'),'file')
     manfile = fullfile(existing_dir,'manifest.txt');
%     fid = fopen(manfile,'r');
%     txt = fread(fid,'uchar=>char')';
%     fclose(fid);
%     infiles = regexp(txt,'\n([\w_.\-]*)\t*0\t*INPUT\t*([\w-]*)','tokens');
%     infiles = cat(1,infiles{:});
%      outfiles = regexp(txt,'\n([\w_.\-]*)\t*0\t*OUTPUT\t*([\w-]*)','tokens');
%      outfiles = cat(1,outfiles{:});
      try
        warning('off','MATLAB:load:variableNotFound')
         for kch = 1:length(files)
            
%             wh= whos('-file',files{kch});
            ldd = load(files{kch},'chan','code'); 
            if ~isfield(ldd,'chan')
                continue
            end
            ldd.chan=ldd.chan(1);
            if ~isfield(ldd.chan,'code')
                ldd.chan.code = '';
            end
          
            ldchin(kch)=ldd;
            
         end
         warning('on','MATLAB:load:variableNotFound')

         chan = [ldchin.chan];
%          if isstruct(chan)
             chann = [chan.channel];
%          end
          [~,unqi]= unique(chann,'stable');
          availch = chan(unqi);
      catch
          floc = fullfile(existing_dir,'model.mat');
           if ~exist(floc,'file')
                    floc = fullfile(existing_dir,'../assets','model.mat');
           end
            ldmod = load(floc);
            availch = ldmod.opts.block.lozchannels;
            flnum = regexp(files,'LFPx(\d*)','tokens','once');
            flnum = cellfun(@str2num,[flnum{:}]);
            files = files(ismember(flnum,[availch.channel]));
                warning('on','MATLAB:load:variableNotFound')

      end
   
%     co2ch([lddat.block.lozchannels.contact]) = [lddat.block.lozchannels.channel];
%     ch2co([lddat.block.lozchannels.channel]) = [lddat.block.lozchannels.contact];
    
    ddat = dir(fullfile(existing_dir,'*hos.mat'));
    dfig = dir(fullfile(existing_dir,'figs'));
    if exist('xne','var') && isa(xne,'xargon') && exist(xne.local_save_dir,'dir') && ~isequal(xne.local_save_dir,existing_dir)
        ddat = [ddat;dir(fullfile(xne.local_save_dir,'*hos.mat'))];
        dfig = [dfig;dir(fullfile(xne.local_save_dir,'figs'))];
    end
    re = regexp({ddat.name},'_(\d*)_hos.*mat','tokens','once');
    if isempty(re)
        doneco = [];
    else
        doneco = cellfun(@str2double,[re{:}]);
    end    
    re = regexp({dfig.name},'contact_(\d*)_bispect.*[.]pdf','tokens','once');
    if ~isempty(re)
        figco = cellfun(@str2double,[re{:}]);
    else
        figco=[];
    end
    missing = ~ismember([availch.contact],doneco);
    
    if ~all(missing)
        lddat = load(fullfile(ddat(1).folder,ddat(1).name),'block','opts','COM');
        if ~isempty(strfind(lddat.COM,'.mat'',ld.blkdat.block,ld.chan.channel);'))
%              ch2co([lddat.block.lozchannels.channel]) = [lddat.block.lozchannels.contact];
%              doneco = ch2co(doneco);
             missing = ~ismember([availch.channel],doneco);
        end

        if ~all(missing) && lddat.opts.make_plots && (~isfield(lddat.opts,'do_regression')||lddat.opts.do_regression)
%             missing = missing | arrayfun(@(x)sum([x.contact]==figco),availch)<lddat.opts.ncomp;
            missing = missing | arrayfun(@(x)sum([x.contact]==figco),availch)<1;
        end
    end
%     if ~all(missing)
% 
%         missing_files = cellfun(@(x)~exist(fullfile(existing_dir,x),'file')&~exist(fullfile(existing_dir,'figs',x),'file'),outfiles(:,1));
%         if exist('xne','var') && isa(xne,'xargon') && exist(xne.local_save_dir,'dir')
%             missing_files = missing_files & cellfun(@(x)~exist(fullfile(xne.local_save_dir,x),'file')&~exist(fullfile(xne.local_save_dir,'figs',x),'file'),outfiles(:,1));        
%         end
% 
%         outfiles(missing_files,:)={'xxxx'};
%     %     donefiles = infiles(ismember(infiles(:,2),outfiles(:,2)),1);
%         if ~isempty(infiles)||isempty(outfiles)
%             [ism,ismi]= ismember(infiles(:,2),outfiles(:,2));
%             if any(ism)
%                 ddat = dir(fullfile(existing_dir,'*hos.mat'));
%                 lddat = load(fullfile(existing_dir,ddat(1).name));
%                 donefiles = infiles(ism);
%                 pdffiles = outfiles(contains(outfiles(:,1),'pdf'));
%                 pdfco = regexp(pdffiles,'contact_(\d*)_bispectral','tokens','once');
%                 pdfco = cellfun(@str2num,[pdfco{:}]);
%                 co2ch([lddat.block.lozchannels.contact]) = [lddat.block.lozchannels.channel];
%                 ch2co([lddat.block.lozchannels.channel]) = [lddat.block.lozchannels.contact];
%     %             [~,c2ci] = ismember(pdfco,[lddat.block.lozchannels.contact]);
%     %             pdfch = [lddat.block.lozchannels(c2ci(c2ci~=0)).channel];
%                 [~,ff,ext] = cellfun(@fileparts,files,'uniformoutput',false);
%                 [fism,ford] = ismember(strcat(ff,ext),infiles(ism,1));
%                 doneco = regexp(outfiles(:,1),'_(\d*)_hos','tokens','once');
%                  [doneco,~,unqi] = unique(cellfun(@str2double,[doneco{ismi(ism)}]),'stable');
%         %         donech = cellfun(@str2double,[donech{ismi(ism)}]);
%             %     files = files(~ismember(strcat(ff,ext),donefiles));
%                 missing = ~ismember(strcat(ff,ext),donefiles);
%             else
%                 missing = true(size(files));
%             end
%         else
%             missing = true(size(files));
%         end       
else
     missing = true(size(files));
     manfile = '';
  
end
% if any(missing)
    xne.jobindices=find( missing );
% else
%    doneco(~missing) =doneco(unqi(ford(~missing)));
% 
%     xne.jobindices = find(~isempty(pdfco) & ~ismember(doneco,pdfco));
% end

%    donef = dir(fullfile(existing_dir,'*_hos.mat'));
%    donech = regexp({donef.name},'_(\d*)_hos[.]mat','tokens','once');
%    donech = cellfun(@str2double,donech);
% %    chs = [opts.block.lozchannels,opts.block.hizchannels];
%    [~,fns] = cellfun(@fileparts,files,'uniformoutput',false);
%    availch = regexp(fns,'(\d*)$','tokens','once');
%    availch = cellfun(@str2double,[availch{:}]);
%    undoneidx = find(~ismember(availch,donech));
blkname = opts.block.block;
summaryfile = [blkname,'_',regexprep(regexp(xne.local_save_dir,['[^',filesep,']*$'],'match','once'),blkname,'')];
xne.finish = @(varargin)finish(xne,summaryfile,varargin{:});

if isempty(xne.jobindices)
   fprintf('\nAll channels finished')
   if ~strcmp(xne.status,'none')
       xne.finish();
   end
   return
elseif xne.nparallel <length(files)
   fprintf('\n%i channels remaining of %i in %s',xne.nparallel,length(files),existing_dir);
end
% end
%     xne.jobindices=undoneidx;  
% else
%     manfile = '';
% end
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

xne.scandep;
xne.dependencies{end+1}=which(opts.stats_function);
xne.dependencies{end+1}=which(opts.plot_function);
deps = [xne.dependencies,matlab.codetools.requiredFilesAndProducts(opts.stats_function)];
deps = [deps,matlab.codetools.requiredFilesAndProducts(opts.plot_function)];
xne.dependencies = unique(deps);
xne.create_job;
% copyfile(which(opts.plot_function),xne.subpaths.mfiles.local)

fid = fopen(fullfile(xne.tempdir,'dependency.m'),'w');
fprintf(fid,'\n%s;',opts.stats_function);
fprintf(fid,'\n%s;',opts.plot_function);
xne.dependencies{end+1} =fullfile(xne.tempdir,'dependency.m'); 

xne.make_bash_script;
xne.make_matlab_wrapper;

if ~isempty(manfile)
    try
        copyfile(manfile,xne.subpaths.output.local);
    end
end



function finish(xne,summaryfile,varargin)

optsfile = fullfile(xne.subpaths.assets.local,'model.mat');
if exist(optsfile)
    if ~exist(xne.local_save_dir,'dir')
        mkdir(xne.local_save_dir)
    end
    copyfile(optsfile,xne.local_save_dir);
end

xne.default_finish();

d = dir(fullfile(xne.local_save_dir,'figs','*contact*cmp*.pdf'));
manfile = fullfile(xne.local_save_dir,'manifest.txt');
if ~isempty(d)
      re = regexp({d.name},'contact_(\d*)_.*cmp(\d*)[.]pdf','tokens','once');
      re = cat(1,re{:});
      cc = cellfun(@str2double,re);
      [srt,srti] = sortrows(cc);
      fns = fullfile(xne.local_save_dir,'figs',{d(srti).name});
elseif exist(manfile,'file')
    fid = fopen(manfile,'r');
    txt = fread(fid,'uchar=>char')';
    fclose(fid);
    re = regexp(txt,'([^\n\s]*[.]pdf)[^\n]','tokens');
    re =[re{:}];
    re = unique(re);
    if isempty(re)
        fprintf('\nNo output files found for %s',xne.local_save_dir);
        return
    end
    re2 = regexp(re,'contact_(\d*)_','tokens','once');
    cnum = cellfun(@(x)str2num(['0',x{:}]),re2);
    [srt,srti] = sort(cnum);
        
    fns = fullfile(xne.local_save_dir,'figs',re(srti(srt>0)));
else
    fns={};
end
if ~isempty(fns)
    com = sprintf('gs -sDEVICE=pdfwrite -dEPSCrop  -dMaxInlineImageSize=100000 -o%s%s%s_summary.pdf %s',xne.local_save_dir,filesep,summaryfile,sprintf(' %s ',fns{:}));
    [err,out]=system(com);
 
    fprintf('%s',out)
    if err~=0
        fprintf('\n\nThe ghostscript command to compile figures into a pdf appears to have failed.\nSee above for clues.\nThe command was as follows:\n\n\t%s',com);
    end
     
end


     
     