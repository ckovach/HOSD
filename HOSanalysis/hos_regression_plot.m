function res=hos_regression_plot(bsidin,useclust,outputdir,outcode)

% Make some plots summarizing the regression result

if nargin < 2
    useclust = false;
end

write_to_wiki=false;

t0=tic;

if nargin < 4 || isempty(outcode)
    outcode = char(java.util.UUID.randomUUID);
end

% if false
%     r = result; %% Here to force dependence
%     h = fact2reg(0);
% end

% addpath signal_toolbox_functions/
%%

% 
% if exist('outputfile','var')
%     ld = load(outputfile);
%     regopts = load(fullfile(inputdir,'regopts'));
%     outputdir2 = fullfile(outputdir,regopts.subdir);
%     if ~exist(outputdir2,'dir')
%         mkdir(outputdir2)
%     end
% else
%     ld= load(fullfile(inputdir,inputfiles{jobindex}));
%     regopts = load(fullfile(inputdir,'regopts'));
%        outputdir2 = fullfile(outputdir,regopts.subdir);
%     if ~exist(outputdir2,'dir')
%         mkdir(outputdir2)
%     end
% %     outputdir2 = outputdir;
% end
% r = regopts.result;

% r.clearfields;

blk = bsidin.block;
% r.subject = blk.subject;
% r.block = blk.block;

% nwin = round(bsidin.fs(1)*regopts.smwint);
% smwin = zeros(size(bsidin.dat));
% h=hann(nwin);
% h=h./sum(h);
% smwin(1:nwin)=h;
% smwin = fft(circshift(smwin,[-floor(nwin/2) 0 ]));
% 
% [T,tt] = chopper(regopts.Trange,regopts.evnt.times,bsidin.fs(1)); 
% T(T<1) = 1;
% T(T>length(ld.dat))=length(ld.dat);
% sz=[length(tt) length(blk.lozchannels)];
    
% hos_permutation_test


% r.contacts = ld.chan.contact;
% r.response_type = regopts.response_label;
% if isfield(ld,'bsids')
%      bsids = ld.bsid;
% else
%     clear bsids
%     xrecs = ld.hos.xrec(ld.dat);
%     ximps = ld.hos.ximp(ld.dat);
%     xfilts = ld.hos.xfilt(ld.dat);
%     
%     stepn = round(ld.hos(1).poverlap*ld.hos(1).buffersize);
%     nget = length(ld.dat) - ld.hos(1).buffersize+1;
%     tindx = (0:ld.hos(1).buffersize-1)';
%     segment(bsdi).wint = (1:stepn:nget);
%    
%     
%      ld.segment(bsdi).tt = fftshift(ld.hos(1).sampt);
% %     segment(bsdi).fs = 1;
% %     segment(bsdi).Trange = [0 ld.hos(1).buffersize];
%     for k = 1:length(ld.hos)
%         bsids(k).xrec = xrecs(:,k);
%         bsids(k).ximp = find(ximps(:,k));
%         bsids(k).xfilt = xfilts(:,k);
%         ld.segment(bsdi).wintadj = ld.segment(bsdi).wint + ld.hos(k).delay;
%         bsids(k).segment = ld.segment;
%         bsids(k).B = ld.hos(k).bicoh;
%         bsids(k).f = ld.hos(k).feature;
%         bsids(k).wb = ld.hos(k).freqindx.Bfreqs{1};
%     end
% end

xrec = bsidin.hos.xrec(bsidin.dat);
ximp = full(bsidin.hos.ximp(bsidin.dat));

for bsdi = 1:length(bsidin.result)
        
%         bsid = bsids(bsdi);
%         
%         switch regopts.measure_type
%             case {'count','impulse'}
%                 imp = zeros(size(bsid.xrec));
%                 imp(bsid.ximp)=1;
%                 cimp = ifft(fft(imp).*conj(smwin));
%                     ttlstring = 'Average Feature Rate';
%                 ylbl = 'rate (1/s)';
%                 minval=0;
%                 scale=bsidin.fs;
%               
%              case 'weighted_impulse'
%                 imp = zeros(size(bsid.xrec));
%                 imp(bsid.ximp)=bsid.xfilt(bsid.ximp);
%                 cimp = ifft(fft(imp).*conj(smwin));
%                     ttlstring = 'Weighted Avgerage Feature Rate';
%                 ylbl = '';
%                 minval=0;
%            case 'rectify'
%                 imp = bsid.xfilt.*(sign(bsid.xfilt)>0);
%                 cimp = ifft(fft(imp).*conj(smwin));
%                 ylbl = '';
%                 minval=0;
%             case 'power'
%                 imp = bsid.xfilt.^2.*(sign(bsid.xfilt)>0);
%                 cimp = ifft(fft(imp).*conj(smwin));
%                 ylbl = 'power';
%                 minval=0;
%             case {'rmcube','cube'}
%                 imp = ifft(fft(bsid.xfilt.^3).*conj(smwin));
%                 if strcmp(regopts.measure_type,'rmcube')
%                     root = 3;
%                     ttlstring = 'Filtered, root-mean cube smoothed';
%                 else
%                     root = 1;
%                     ttlstring = 'Filtered, cube smoothed';
%                 end                        
%                 cimp = sign(imp).*abs(imp).^(1/root);
% 
%                 ylbl = 'root-mean cube';
%                 minval=-1;
%             otherwise
%                 error('Unrecognized measure type');
%         end
% %         cimpc = cimp;
%         intercept = true;
%         if isfield(regopts,'center') && regopts.center || length(regopts.eReg)==1 && size(regopts.eReg.value,2)==1&& all(abs(1-regopts.eReg.value)<=eps)
%             cimpc = cimp(T)-repmat(mean(cimp(T)),size(T,1),1);
%             intercept = false;
%             glmargs = {'intercept',false};
%         else
%             cimpc = cimp(T);
%             glmargs = {};
%         end
%         [a,b,c,iXX,serr] = complexglm(cimpc',[regopts.eReg.value],glmargs{:});
% %         discard_times = std(a)*serr<eps;  %%% The gaussian model is inapproprite in these regions, so they will be discardexd.
%           se= sqrt(diag(iXX)*serr);
%          wald = a./se;% + 0./~discard_times;
% %          b = b;% + 0./~discard_times;
%          
%          no_imp_trial =std(cimp(T)).*std(cimp)<=eps;
%          if size(wald,1)>1 && intercept
%              wald(end,:)=[]; % Ignore baseline
%          end
% %             if length(regopts.eReg)>1
% %                 [a0,b0,c0] = complexglm(cimp(T)',[regopts.eReg(end).value]);
% %             else
%             b0=0;
% %             end
%         if exist('permtest','var')
%             ptrans = permtest.bsp_permP(:,bsdi); 
%         
%         else
%             ptrans = 1-chi2cdf(b-b0,sum([regopts.eReg.Npar]));
%         end
%          if ~all(isnan(ptrans))
%             qtrans = fdr(ptrans);
%         else
%             qtrans = ptrans;
%          end
%             
%         res.Pbsp(:,bsdi)=ptrans;
%         res.Qbsp(:,bsdi)=qtrans;
        xr = xrec(:,bsdi);
        xi = ximp(:,bsdi);
        fig = figure;
        
        dbx = dbt(bsidin.dat,bsidin.fs,4,'upsample',2);

        if isfield(bsidin.result(bsdi),'fit') && ~isempty(bsidin.result(bsdi).fit)
            mdl=bsidin.result(bsdi).fit.model;
            evw = mdl.get_event_window;
            tt = evw.tt;
            T = evw.T;
            [unqev,~,unqevi] = unique(mdl.event.evnt','rows');
            [srt,srti] = sort(unqevi);
           if isnumeric(unqev)
                unqev = arrayfun(@(k)sprintf('%i ',unqev(k,:)),1:size(unqev,1),'uniformoutput',false);
            end
            fit = bsidin.result(bsdi).fit;
            windest = [fit.regressors.windowest];
    %         smwin = hann(round(.1*bsidin.fs));
    %         smwin = smwin./sum(smwin);
            xism = convn(xi,ones(round(.05*bsidin.fs)),'same')>0;
    %         xism2 = convn(xi,smwin,'same');
            for evk = 1:length(unqev)
                      res.Mevk(:,evk,bsdi) = mean(xism(T(:,evk==unqevi)),2);
    %                  res.Mevk(:,evk,bsdi) = fit.regressors;

                    res.Merp(:,evk,1) = mean(bsidin.dat(T(:,evk==unqevi)),2);

                    res.compMerp(:,evk,bsdi) = mean(xr(T(:,evk==unqevi)),2); %component contribution to the erp

            end
               windowreg=  find(arrayfun(@(x)~isempty(x.windowest),fit.regressors));
            ax1 = subplot(2+length(windowreg),3,1);

            imh = imagesc(tt,[],xism(T(:,srti))');
            edges = [1;find(diff(srt)>0);size(T,2)];
            hold on, plot(tt([1 end]),[1 1]'*edges','w--');
            set(ax1,'ytick',(edges(1:end-1)+edges(2:end))/2,'yticklabel',unqev);
            title('Raster')

            ax1 = subplot(2+length(windowreg),3,4);

            plh = plot(tt,res.Mevk(:,:,bsdi));
            title('Mean Raster')
    %         legend(unqev)
            axis tight
            grid on
            for k = 1:length(windowreg)
                subplot(2+length(windowreg),3,3*k+4)
                rg = fit.regressors(windowreg(k));
                plot(tt([1 end]),[0 0],'k')
                hold on
                plot(tt,rg.windowest.wald);
                axis tight
                ylim([-1 1]*max([abs(ylim) 4]))
                tlt = sprintf('%s: ',rg.label);
                if rg.waldpval<.001
                    tlt = sprintf('%s*** P= %0.2g',tlt,rg.waldpval);
                elseif rg.waldpval<.001
                    tlt = sprintf('%s** P= %0.2g',tlt,rg.waldpval);
                elseif rg.waldpval<.05
                    tlt = sprintf('%s* P= %0.2g',tlt,rg.waldpval);
                else
                    tlt = sprintf('%sn.s. P= %0.2g',tlt,rg.waldpval);
                end

                title(tlt)
                grid on
    %             yl=ylim;
    %             set(gca,'ytick',floor(min(yl)):ceil(max(yl)))

            end
            xlabel('s')
            ylabel('est./stderr.')
            
            
        dm  = mdl.designMtx;
        X = [dm.value];
     
        kp = arrayfun(@(x)x.levmat(end,:)==1&~isempty(x.window),dm,'uniformoutput',false);
        X = X(T(round(end/2),:),[kp{:}]);
        if ~isfield(mdl.event,'center')||mdl.event(1).center
            X(:,end)= [];
        end
%         X = X(T(round(end/2),:),all(dm(1).levmat(end,:)==1,1));
         [a2,b2,c2] = complexglm(bsidin.dat(T)',X);
%              if length(regopts.eReg)>1
%                 [a20,b20,c20] = complexglm(ld.dat(T)',[regopts.eReg(end).value]);
%              else
             b20=0;
%              end
        if exist('permtest','var')
            ptrans2 = permtest.erp_permP; 
        
        else
             ptrans2 = 1-chi2cdf(b2-b20,size(X,2));
        end
        if ~all(isnan(ptrans2(:)))
            qtrans2 = ptrans2;
            qtrans2(:) = fdr(ptrans2(:));
        else 
            qtrans2 = ptrans2;
        end
        
        ax2 = subplot(4,3,8);
        
        plh2 = plot(tt,squeeze(res.Merp));
        hold on
        plh3 = plot(tt,squeeze(res.compMerp(:,:,bsdi)),'linewidth',2);
%          ax2.Position(3)=.5;
        res.Perp(:,bsdi,:)=ptrans2;
          res.Qerp(:,bsdi,:)=qtrans2;
%            if length(regopts.eReg)==2
%                set(plh2(2:2:end),'color',[1 1 1]/2);
%                  set(plh2(1:2:end),'color','r');
%            else
%                for plk = 1:length(plh2)
%                    set(plh2(plk),'color',cols(plk,:))
%                end
%            end
        
       axis tight
           hold on, plot(tt,0./(qtrans2<1),'k','linewidth',.5)
%         plot(tt ,0./(qtrans2<.1),'k','linewidth',2)
        plot(tt ,0./(qtrans2<.05),'k','linewidth',2)
        plot(tt  ,0./(qtrans2<.01),'k','linewidth',4)
        hold on,th= plot([0 0],ylim,'k');
%         th(2)=plot([1 1],ylim,'k--');
%         th(3)=plot([2 2],ylim,'k:');
        grid on
        [A2,att2] = choptf(mdl.event(1).Trange,mdl.event(1).times(:)',dbx,mdl.event(1).Trange);       
        AA  = reshape(20*log10(abs(A2)),numel(A2(:,:,1)),size(A2,3));
        [a3,b3,c3] = complexglm(AA',X,'diagonly',false);
%             if length(regopts.eReg)>1
%             [a30,b30,c30] = complexglm(AA',[regopts.eReg(end).value]);
%             else
            b30=0;
%             end
        ptrans3 = 1-chi2cdf(b3-b30,size(X,2));
        qtrans3 = reshape(ptrans3,size(A2(:,:,1)));
        if ~all(isnan(ptrans3(:)))
            qtrans3(:) = fdr(ptrans3(:));
        end
        res.Perbp(:,bsdi)=ptrans3;
            res.Qerbp(:,:,bsdi)=qtrans3;
            res.att2=att2;
            res.freq=dbx.frequency;
%             if length(regopts.eReg)==2
%     %         lg= legend([plh2(1:2);th'],{sprf(regopts.unqev(1:2:end)),sprf(regopts.unqev(2:2:end)),'Onset','Transition','End'},'position',[ 0.0053    0.5696    0.1652    0.1237]);
         lg= legend(plh2,unqev,'position',[ 0.69   0.75   0.15   0.1]);
%              cols = jet(length(plh2)).*.9;
%         cols = hsv2rgb([(0:length(regopts.unqev)-1)'/4,ones(4,2)*.9]);

         for lgi = 1:length(plh2)
%              plh2(lgi).Color=cols(lgi,:);
              plh3(lgi).Color=plh2(lgi).Color;
         end
%              lg= legend([plh2;th'],{regopts.eReg.label,'Onset'},'position',[ 0.0053    0.5696    0.1652    0.1237]);
%             elseif length(regopts.eReg)==3
%                   mkrow=@(x)x(:)';
%                   lg= legend([plh2([2,1:2:end]);th'],[{sprf(control_events)},arrayfun(@num2str,mkrow(regopts.unqev(mod(regopts.unqev,10)==0)),'uniformoutput',false),{'Onset'}],'position',[ 0.0053    0.5696    0.1652    0.1237]);
% 
%             else
%                 %         lg= legend([plh2;th'],[arrayfun(@num2str,regopts.unqev(:)','uniformoutput',false),{'Onset','Transition','End'}],'position',[ 0.0053    0.5696    0.1652    0.1237]);
%               mkrow=@(x)x(:)';
%             lg= legend([plh2;th'],[{sprf(control_events)},arrayfun(@num2str,mkrow(regopts.unqev(regopts.unqev>1)),'uniformoutput',false),{'Onset'}],'position',[ 0.0053    0.5696    0.1652    0.1237]);
%             end
            title('Stim. evoked response & component')

          ax5=subplot(4,3,11);
          getfr = find(dbx.frequency>0);
          wintp =10.^(linspace(log10(dbx.frequency(getfr(1))),log10(dbx.frequency(end)),length(dbx.frequency)));   
          MA = mean(20*log(abs(A2)),3)';
          Aintrp = interp2(att2,dbx.frequency(getfr),MA(getfr,:),att2,wintp);    
          imagesc(att2,log10(wintp),Aintrp)
              set(ax5,'ytick',log10(2.^(0:log2(wintp(end)))),'yticklabel',2.^(0:log2(wintp(end)))) %     caxis([0 quantile(BB(:),.999)])

    %       imagesc(att2,dbx.frequency,mean(20*log(abs(A2)),3)')
          caxis([-1 1]*min(abs(caxis)))
              hold on
            [~,ch] =contour(att2,log10(dbx.frequency(getfr)),-log(qtrans3(:,getfr)'),-log([.05 .05])); 
    %         [~,ch(2)] =contour(att2,dbx.frequency,-log(qtrans3'),-log([ .05 ])); 
    %         [~,ch(3)] =contour(att2,dbx.frequency,-log(qtrans3'),-log([.01])); 
            set(ch,'color','black','linewidth',1)
    %         set(ch(1),'linestyle',':')
    %         set(ch(2),'linewidth',1)
    %         set(ch(3),'linewidth',2)


          axis xy
            title('Stim. induced')
    %            colorbar('north')
    %       colormap jet
          set(fig,'position',[ 584    30   696   671]);
          xlabel('trial time (s)')
            ylabel Hz
        end
%         hold on, plot(tt,0./(qtrans<1),'k','linewidth',.5)
% 
%         plot(tt ,0./(qtrans<.05),'k','linewidth',2)
%         plot(tt  ,0./(qtrans<.01),'k','linewidth',4)
%         ylim(1.1*[minval 1]*max(ylim))
% 
%         hold on, th=plot([0 0],ylim,'k');

%         grid on
%          ax1.Position(3)=.5;

%         ax= axis;

%         ttl=title(ttlstring);

%            ylabel(ylbl)
%             ylabel('$\frac{\hat{\beta}}{\mathrm{std.\;err.}}$','fontsize',20,'interpreter','latex')
%                 cols = jet(length(plh));
%                 cols = hsv2rgb([(0:length(regopts.unqev)-1)'/4,ones(4,2)*.9]);
%             cols = hsv2rgb([(0:length(regopts.unqev)-1)'/length(regopts.unqev),ones(length(regopts.unqev),2)*.9]);
%          for lgi = 1:length(plh)
%              plh(lgi).Color=cols(lgi,:);
%          end
%             xlabel('time (s)')
%         subi = @(x)x(1:end-1);
%         sprf =@(x)subi(sprintf('%i,',x));
%         legend([plh(1:2);th'],{sprf(regopts.eventids(1:2:end)),sprf(regopts.eventids(2:2:end)),'Onset','Transition','End'},'location','northwest')


% %%%%%%%%%%%%%% Effect plot
%         ax3 =subplot(5,3,5);
% %             plh = plot(tt,squeeze(res.Mevk(:,kk,:,bsdi)));
%         plh = plot(tt,squeeze(wald));
% 
%         hold on, plot(tt,0./(qtrans<1),'k','linewidth',.5)
%         plot(tt ,0./(qtrans<.05),'k','linewidth',2)
%         plot(tt  ,0./(qtrans<.01),'k','linewidth',4)
%           ylim(1.1*[-1 1]* max(max(abs(wald(:))),4))
%         hold on, th=plot([0 0],ylim,'k');
%         grid on
%      ax3.Position(3)=.5;
%         ax= axis;
%         ttl=title('Effects (Wald statistic)');
%         ylabel('$\frac{\hat{\beta}}{\mathrm{std.\;err.}}$','fontsize',16,'interpreter','latex')
% 
%         xlabel('time (s)')
%         subi = @(x)x(1:end-1);
%         sprf =@(x)subi(sprintf('%i,',x));
%         lg= legend([plh],{regopts.eReg.label,'Onset'},'position',[ 0.67    0.454    0.233    0.085]);

%%%%%%%%%%%%%%%%%%

        subplot(4,3,6);
        plot(ifftshift(bsidin.hos(bsdi).sampt)/bsidin.fs,bsidin.hos(bsdi).feature)
        axis tight
        title(sprintf('Component %i waveform',bsdi))
       grid on
        subplot(4,3,9);
        wwin= (bsidin.hos(bsdi).sampt)/bsidin.hos(bsdi).buffersize*bsidin.fs;
        plot(wwin,abs(fft(bsidin.hos(bsdi).xfilt(bsidin.hos(bsdi).feature))))
        xl = xlim;
        axis tight
        xlim([2*bsidin.fs/bsidin.hos(bsdi).buffersize xl(2)] )
        set(gca,'xscale','log','xtick',2.^(0:log2(max(wwin))))
        title(sprintf('Component %i Normalized FFT',bsdi))
       grid on

       [A,att] = choptf(bsidin.segment(bsdi).Trange*bsidin.segment(bsdi).fs/bsidin.fs,bsidin.segment(bsdi).wintadj*bsidin.segment(bsdi).fs/bsidin.fs,dbx,bsidin.segment(bsdi).Trange*bsidin.segment(bsdi).fs/bsidin.fs);
        ax6=subplot(4,3,12);
       getfr = find(dbx.frequency>0);
       wintp =10.^(linspace(log10(dbx.frequency(getfr(1))),log10(dbx.frequency(end)),length(dbx.frequency)));   
       MA = mean(20*log(abs(A)),3)';
       Aintrp = interp2(att,dbx.frequency(getfr),MA(getfr,:),att,wintp);    
    
%        imagesc(att,dbx.frequency,mean(20*log10(abs(A)),3)')
        imagesc(att,log10(wintp),Aintrp)
        set(ax6,'ytick',log10(2.^(0:log2(wintp(end)))),'yticklabel',2.^(0:log2(wintp(end)))) %     caxis([0 quantile(BB(:),.999)])
         axis xy
        title('Feat. induced ')
         caxis([-1 1]*min(abs(caxis)))
         xlabel('(s)')
         
         
        
          
        
       wb = fftshift(bsidin.hos(1).freqindx.Bfreqs{1});
        BB=fftshift(abs(bsidin.hos(bsdi).bicoh));
        BB(wb<=0,:)=[];
        BB(:,wb<=0)=[];
        PB=fftshift(abs(bsidin.hos(bsdi).partialbicoh));
        PB(wb<=0,:)=[];
        PB(:,wb<=0)=[];
        
        wb(wb<=0)=[];
         wintp =10.^(linspace(log10(wb(2)),log10(wb(end)),length(wb)));
        mm = meshgrid(wb,wb);
        mmintp = meshgrid(wintp,wintp);
        BBintp = interp2(mm,mm',BB,mmintp,mmintp');
        PBintp = interp2(mm,mm',PB,mmintp,mmintp');
        
        BBintp(mm<mm')=PBintp(mm<mm');
        ax4 = subplot(2,3,2);

%        imagesc(wb,wb,abs(ld.bsid(1).B))
       imagesc(log10(wintp),log10(wintp),BBintp)
          set(ax4,'xtick',log10(2.^(0:log2(wintp(end)))),'xticklabel',2.^(0:log2(wintp(end))),...
                  'ytick',log10(2.^(0:log2(wintp(end)))),'yticklabel',2.^(0:log2(wintp(end)))) %     caxis([0 quantile(BB(:),.999)])
       axis image xy
       cbar = colorbar('SouthOutside','position',[  0.4303    0.5488    0.1767    0.0200]);
%           axis([0 1 0 1]*max(wb))
%          ax4.Position = [.71 .53 .188 .4];
      title('Part/Full Bicoherence');
      xlabel Hz
      ylabel Hz
%       if bsdi==1
          cax = caxis;
%       else
%           caxis(cax);
%       end
      set(ax4,'xtick',get(ax4,'ytick'))
      
%         ax4 = subplot(4,3,5);
%        imagesc(log10(wintp),log10(wintp),PBintp)
%           set(ax4,'xtick',log10(2.^(0:log2(wintp(end)))),'xticklabel',2.^(0:log2(wintp(end))),...
%                   'ytick',log10(2.^(0:log2(wintp(end)))),'yticklabel',2.^(0:log2(wintp(end)))) %     caxis([0 quantile(BB(:),.999)])
%        axis image xy
%        cbar = colorbar;
%        caxis(cax);
% %           axis([0 1 0 1]*max(wb))
% %          ax4.Position = [.71 .53 .188 .4];
%       title('Partial Bicoherence');
%       xlabel Hz
%       ylabel Hz
%       set(ax4,'xtick',get(ax4,'ytick'))
% %           caxis([0 1])

      ax0 = axes;axis off
       tth = text(.728,.975,sprintf('Contact %i:%s %i\nComponent: %i',bsidin.chan.contact,regexprep(bsidin.chan.label,'_','\\_'),bsidin.chan.number,bsdi));
%             tth.Positin = [-2 .05
       %        tth = title(sprintf('Contact %i: %s %i',ld.chan.contact,bsidin.chan.label,bsidin.chan.number));
       tth(2) = text(.728,1.051,sprintf('Block: %s',bsidin.block.block));
       tth(3) = text(.728,1.02,sprintf('Protocol: %s',bsidin.block.subprotocol));
%            axis off

%       c = r.available_contacts([r.available_contacts.Contact_Number]==r.contacts);

%       hl = strcat(sprintf('%sindex.php/',r.url),{c.main,bsidin.block.block,bsidin.block.subprotocol});
        hl = '';
        if useclust
%          locfn= fullfile(bsidin.opts.outputdir2,sprintf('temp%icmp%i.pdf',r.contacts,bsdi));
%          pdflink(fig,tth,hl,locfn)
          figdir = fullfile(outputdir,'figs');
         if ~exist(figdir,'dir')
             mkdir(figdir)
         end
         wfnpart = sprintf('%s_contact_%03i_bispectral_%s_cmp%i.pdf',bsidin.block.block,bsidin.chan.contact,regexprep(bsidin.block.subprotocol,'\s','_'),bsdi);
         wfn = fullfile(figdir,wfnpart);
          pdflink(fig,[],[],wfn)
                    
          fid = fopen(fullfile(outputdir,'manifest.txt'),'a+');
          fprintf(fid,'\n%s\t0\tOUTPUT\t%s\t%0.3fs',wfnpart,outcode,toc(t0));
          fclose(fid);

%       pdflink(fig,tth,hl,wfn)

         %          descr = sprintf('Comparision of bispectral feature analysis and ERP for the figure ground task: block %s, contact %i',bsid.block.block,r.contacts);
%         r.pdf=wfn;
%         url = {regopts.url};%,regopts.url2};
%         rev = {regopts.rev,regopts.rev};
%         cc = cat(1,url,rev);
%         str = sprintf(',%s',url{:});
%         str(1)='';
%         r.repository_url=str;
% %         r.revision_number=rev{1};
%          r.analysis='Bispectral feature';

%              se= sqrt(diag(iXX)*serr);
%              wald = a./se;
%          [mn,mni] = min(ptrans' +0./(tt>=-Inf & tt <=Inf));
%          [~,mxiw] = max(abs(wald(:,mni)))
%              r.values = wald(1:2,mni);
%              r.values = -log10(mn);
%          r.pvalue = min(qtrans);
%          if min(qtrans)<.05
%              eff=' ''''''significant''''''';
%          else
%               eff=' non-significant';
% 
%          end
% %          if wald(1,mni)<0
% %             effect = [eff,' decrease'];
% %          else
% %              effect = [eff,' increase'];
% % 
% %          end
%               descr = sprintf('=\n===%s: Contact %i, Component %i===\nThe comparison for component %i was associated with %s modulation of the feature returned by bispectral feature extraction (minimum FDR corrected value P=%0.3g at t= %gs).\nMin q for ERP:%0.2g\nMin q for ERBP:%0.2g.',...
%                   bsidin.block.subject,bsidin.chan.contact,bsdi([1 1]),eff,qtrans(mni),tt(mni),min(qtrans2(:)),min(qtrans3(:)));

% 
%          r.description = descr;
%          r.pdf = wfn;
%           r.page_name = sprintf('%s/Contact %i %s Result Component %i',bsidin.block.block,bsidin.chan.contact,bsidin.block.subprotocol,bsdi);
%           res.upload_error=false;
%           if write_to_wiki && min([qtrans(:);qtrans2(:);qtrans3(:)])<regopts.qthresh
%             try
%               r.upload(locfn,[],'pdf',descr,wfn,regopts.clobber);
%               r.write(regopts.clobber);
%             catch err
%                 res.upload_error = err;
%             end
%           elseif write_to_wiki
%                 txt = r.readPage(r.page_name);
%              if ~isempty(txt)
%                 r.query_params=struct('action','delete','title',r.page_name,'reason','Deleted after non-significant effect in the updated analysis','token',r.get_token);
%                 q=r.request;
%              end
%           end
          %%
          res.pdfpages{bsdi}=wfn;

        end
end
% res.feat = [bsids.f];
fid = fopen(which(mfilename),'r');
res.COM = fread(fid,'uchar=>char')';
% res.opts = regopts;
res.block=bsidin.block;
res.chan=bsidin.chan;
fclose(fid);

 
% if useclust
%     save(fullfile(outputdir,sprintf('%s_%s_%s_ch%i%s',blk.block,regopts.label,mfilename,bsidin.chan.channel,regopts.sufx)),'-struct','res')
% end


%     if err==0
%         r.upload(outfn,[],'pdf',sprintf('Summary of results for response [[Summary of::%s]] for block [[Has block::%s]]',r.response_type,blk.block),regopts.clobber);
%         system(sprintf('rm -r %s',tmpdir))
%     end
    
    