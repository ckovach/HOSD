function varargout =plotpolycoh(hos,ax,plotwhat,absfun)

%Plot polycoherence. For orders greater than 3, plot 2d slices.


if nargin < 2 
    ax = [];
end
thresh = 0;
inputfun = @(x)x;
if nargin < 3 || isempty(plotwhat)
    plotwhat = 'abs';
elseif isnumeric(plotwhat) && isscalar(plotwhat)
    thresh = plotwhat;
    plotwhat = 'phase';
else
    inputfun = @(x)plotwhat(x);
    plotwhat = 'input';
end
%scale = 'log';
%   scale = 'lin';
 scale = '';
switch scale

    case 'log'
        wtr = @(x)10.^x;
        iwtr = @(x)log10(x);
    case {'lin','linear'}
        wtr = @(x)x;
        iwtr = wtr;
    otherwise
        wtr = @(x)x;
        iwtr = wtr;
end

plot_partialbc = false;

% fig = figure('WindowButtonMotionFcn',@(a,b,c)figcallback(a,wtr,iwtr));
% fig = figure;
switch plotwhat
    case {'mag','abs','magnitude'}
        absfun = @abs;
    case {'phase','angle'}
      absfun = @(x)angle(x) + 0./(abs(x)>thresh);
    case 'input'
        if nargin < 4
          absfun = @(x)x;
        end
    otherwise
        if isa(plotwhat,'function_handle')
            if nargin < 4
                absfun = plotwhat;
            end
                
        else
            error('Unrecognized plotting option')
        end
end
%absfun = @real;
subxy = almostSquare(length(hos));
txth = [];
for hi = 1:length(hos)
   
    
    if hos(hi).diagonal_slice && hos(hi).order == 4
        subxy = [length(hos) 2];
        wb0 = cellfun(@fftshift,hos(hi).freqindx.Bfreqs,'uniformoutput',false);
     
        power = wb0{1};
        power = power(power>0 & power >hos(hi).highpass &power< hos(hi).lowpass);
        modfreq = wb0{end}-min(wb0{end});
        modfreq = modfreq(abs(modfreq)>hos(hi).shighpass & abs(modfreq)<=hos(hi).slowpass);
        
        [P,M] = meshgrid(power,modfreq);
        [W1,W2,W3] = meshgrid(wb0{:});
        B = nan*M;
        Bpart = B;
        B(:) = interp3(W1,W2,W3,fftshift(inputfun(hos(hi).bicoh)),P(:),P(:),-P(:)+M(:));
        Bpart(:) = interp3(W1,W2,W3,fftshift(inputfun(hos(hi).partialbicoh)),P(:),P(:),-P(:)+M(:));
        Bpart(isnan(B))=nan;
         if hi>size(ax,1)
             ax(hi,1) = subplot(subxy(1),subxy(2),1 + 2*(hi-1));
             ax(hi,2) = subplot(subxy(1),subxy(2),2 + 2*(hi-1));
%         else
%             axes(ax(hi))
         end
        imh = pcolor(power,modfreq,absfun(B),'parent',ax(hi,1));
        set(imh,'facecolor','flat','edgecolor','none');
        if hi==1
         title(ax(hi,1),sprintf('Modulogram (Trispectrum Diagonal Slice)'))
        else
         title(ax(hi,2),sprintf('Residual Modulogram After Comp. %i',hi-1))
        end    
        axis xy
%         xlabel('Envelope modulation freq.(Hz)')
%         ylabel('Band freq. (Hz)')
        ylabel('Envelope modulation freq.(Hz)')
        xlabel('Band freq. (Hz)')
        
         if hi==1
            cax = caxis;
        else
            caxis(cax)
         end
        imh = pcolor(power,modfreq,absfun(Bpart),'parent',ax(hi,2));
        set(imh,'facecolor','flat','edgecolor','none');
         title(ax(hi,2),sprintf('Component %i Trispectrum Modulogram',hi))
        axis xy

    elseif hos(hi).diagonal_slice && hos(hi).order == 3
        wb0 = cellfun(@fftshift,hos(hi).freqindx.Bfreqs,'uniformoutput',false);
        B = nan*wb0{1};
        Bpart = B;
        [W1,W2] = meshgrid(wb0{:});
        B(:) = interp2(W1,W2,fftshift(absfun(inputfun(hos(hi).bicoh))),wb0{1},wb0{1});
        Bpart(:) = interp2(W1,W2,fftshift(absfun(inputfun(hos(hi).partialbicoh))),wb0{1},wb0{1});
         if hi>length(ax)
             ax(hi) = subplot(subxy(1),subxy(2),hi);
%         else
%             axes(ax(hi))
         end
        imh = plot(wb0{1},B.*(wb0{1}>=0)+ (1-(wb0{1}>=0)).*Bpart,'parent',ax(hi));
         title(sprintf('Component %i Bispectrum Diagonal Slice',hi))
         xlabel('freq.(Hz)')
    else
        wb0 = cellfun(@fftshift,hos(hi).freqindx.Bfreqs,'uniformoutput',false);
        do_interp=true;
        switch scale

            case 'log'
                wb = cellfun(@(x)log10(logspace(log10(hos(hi).highpass),log10(max(x)),length(x))),wb0,'uniformoutput',false);
                tick = @(x)(linspace(x(1),x(end),10));
                ticklbl = @(x)round(10.^tick(x)./10.^(round(tick(x))-1)).*10.^(round(tick(x))-1);
                wbin = cellfun(@(x)log10(abs(x)),wb0,'uniformoutput',false);
  
            case {'lin','linear'}
                wb = cellfun(@(x)x(x>hos(hi).highpass),wb0,'uniformoutput',false);
                tick = @(x)(linspace(x(1),x(end),10));
                ticklbl = @(x)tick(x);
               wbin = wb;
                
            otherwise
                wb=wb0;
                do_interp=false;
                wbin = wb;
  
        end
        [wbin{:}] = ndgrid(wb0{:});
      
        wb(end+1) = {0};

%         ax(hi) = subplot(subxy(1),subxy(2),hi,'UserData',wb);
        if hi>length(ax)
             ax(hi) = subplot(subxy(1),subxy(2),hi,'UserData',wb);
%         else
%             axes(ax(hi))
         end
        hold on

        dwb = cellfun(@(x)diff(x([1 end])),wb(1:2));

        nsamps = cellfun(@length,wb);

        K = cellfun(@(x)1:length(x),wb(3:end),'uniformoutput',false);
        
        Ws = wb;
        [K{:}] = ndgrid(K{:});
        [Ws{:}] = ndgrid(wb{:});
        Ws{end} = -sum(cat(length(Ws),Ws{1:end-1}),length(Ws));
        if length(wb)>3
                K = cellfun(@(x)x(:),K,'uniformoutput',false);
                K = [K{:}];
                K = K(K(:,end-1)<=length(wb{end-1})/2,:);
        else
            K=1
        end
        dims = [];
        for k = 3:length(wb)
    %         rg = [0:sqrt(length(wb{k}))];
    %         divs = arrayfun(@(x)find(rem(x,1:x)==0),nsamps(k)+rg,'uniformoutput',false);
    %         divr = arrayfun(@(x,y)y{1}-x./y{1},nsamps(k)+rg,divs,'uniformoutput',false);
    %         
    %         [mn,mni] = cellfun(@(x,y)min(abs(x)),divr,divs);
    %         [~,mni2] = min(mn);
    %         dm = divs{mni2}(mni(mni2));
            dims(k-2,:) = almostSquare(nsamps(k));
        end    
%         wbin = wb(1:end-1);
        wbout=wb(1:end-1);
        [wbout{:}] = ndgrid(wb{:});

        B = interpn(wbin{:},absfun(fftshift(hos(hi).bicoh)),wbout{:});
%         B = fftshift(full(hos(hi).freqindx.keep));
        if plot_partialbc 
        Bpart = interpn(wbin{:},absfun(fftshift(hos(hi).partialbicoh)),wbout{:},'parent',ax(hi));
        else
            Bpart = B;
        end
        for k = 1:size(K,1)

            cent = [0 0];
            for kk = 1:size(K,2)
                cent =  cent + [mod(K(k,kk)-1,dims(kk,2)), floor((K(k,kk)-1)/dims(kk,2))].*prod(dims(1:kk-1,:),1);
            end        
            wcent = cent.*dwb;

            imk = num2cell(K(k,:));
            Bplot = B(:,:,imk{:}).*(wbout{2}(:,:,imk{:})>=wbout{1}(:,:,imk{:})) + (1-(wbout{2}(:,:,imk{:})>=wbout{1}(:,:,imk{:}))).*Bpart(:,:,imk{:});
            imh(imk{:}) = imagesc(wb{1}+wcent(1),wb{2}+wcent(2),Bplot,'parent',ax(hi));
            delete(datatip(imh(imk{:})));
%             delete(dt)
            for kk = 1:length(Ws)
                imh(imk{:}).DataTipTemplate.DataTipRows(kk) = dataTipTextRow(sprintf('f_%i:',kk),Ws{kk}(:,:,imk{:}));
            end
                imh(imk{:}).DataTipTemplate.DataTipRows(length(Ws)+1) = dataTipTextRow('Value:',B(:,:,imk{:}));
            for kk = 3:length(wb)-1
                if kk==3
                    txt{kk-2} = sprintf('%0.2f',wb{kk}(K(k,kk-2)));
                else
                    txt{kk-2} = sprintf(', %0.2f',wb{kk}(K(k,kk-2)));
                end
            end
            if ~isempty(kk)
             txth(imk{:})=text(wcent(1)+dwb(1)/2+wb{1}(1),wcent(2)+.9*dwb(2)+wb{2}(1),[txt{:}],'Color','w','HorizontalAlignment','center','fontsize',8,'parent',ax(hi));
            end
        end
        title(sprintf('Component %i',hi))
        axis image
        if hi==1
            cax = caxis;
        else
            caxis(cax)
        end
%         tick = @(x)round(linspace(x(1),x(end),10));
%          set(gca,'xtick',tick(wb{1}),'ytick',tick(wb{2}),'xticklabel',ticklbl(wb{1}),'yticklabel',ticklbl(wb{2}))
    end
    if nargout >=3
        Bs{hi} = B;
    end
end

colorbar

if nargout>=1
    varargout{1}=imh;
end
if nargout>=2
    varargout{2} = txth;
end
if nargout >=3
    varargout{3} = Bs;
end

% function figcallback(ax,wtr,iwtr)
%     
% ax = gca;
% wb = ax.UserData;
% lblx = @(x)wtr(mod(x,iwtr(max(wb{1}))));
% lbly = @(x)wtr(mod(x,iwtr(max(wb{2}))));
% 
% set(ax,'xticklabel',lblx(get(ax,'xtick')),'yticklabel',lblx(get(ax,'ytick')));
%     
    

function out =  almostSquare(N)

rg = 0:sqrt(N);
divs = arrayfun(@(x)find(rem(x,1:x)==0),N+rg,'uniformoutput',false);
divr = arrayfun(@(x,y)y{1}-x./y{1},N+rg,divs,'uniformoutput',false);

[mn,mni] = cellfun(@(x,y)min(abs(x)),divr,divs);
[~,mni2] = min(mn);
dm = divs{mni2}(mni(mni2));
out = [(N+rg(mni2))/dm  dm];
