function varargout = modgram(me,Bplot0)

%Extract the modulogram from a 4th order spectrum

if me(1).order~=4
    error('Spectrum is not 4th order');
end

if nargin > 1 && length(Bplot0) == length(me.B)
    Bplot0 = Bplot0(me.fullmap);
end

fun = @(x)x;

for cmpi = 1:length(me)
    if nargin < 2
        B0 = fftshift(me(cmpi).bicoh);
        inds = [repmat({':'},1,me.order-1),{1}];
        Bplot =  B0(inds{:});
        for k = 2:size(B0,me(cmpi).order)
            inds(me.order) = {k};
             Bplot = cat(me(cmpi).order,Bplot,fftshift(B0(inds{:})));
        end
        fun = @abs;
    elseif isa(Bplot0,'function_handle')
        Bplot = fftshift(me(cmpi).bicoh);
        fun = Bplot0;
    elseif size(Bplot0,2) == length(me)
        Bplot = Bplot0(:,cmpi);
        Bplot = fftshift(Bplot(me(cmpi).fullmap));
    else 
        inds = [repmat({':'},1,me.order-1),{1}];
        B0 = Bplot0;
        Bplot =  fftshift(B0(inds{:}));
    end

    wb0 = cellfun(@fftshift,me(cmpi).freqindx.Bfreqs,'uniformoutput',false);
    [W1,W2,W3] = meshgrid(wb0{:});
     
    power = wb0{1};
   % power = power(power>0 & power >me(cmpi).highpass &power< me(cmpi).lowpass);
    power = power(power>0);
    modfreq = wb0{end}-min(wb0{end});
    modfreq = modfreq(abs(modfreq)<me.slowpass);

  %  modfreq = modfreq(abs(modfreq)>me(cmpi).shighpass & abs(modfreq)<=me(cmpi).slowpass);
    modfreq = modfreq(abs(modfreq)>0);
    
    [P,M] = meshgrid(power,modfreq);
    Bout =[];
    for k = 1:size(Bplot,me(cmpi).order)
        B = nan*M;
        Bpart = B;
         inds{me.order} = k;
   
        B(:) = interp3(W1,W2,W3,Bplot(inds{:}),P(:),P(:),-P(:)+M(:),'nearest');
        Bout = cat(me(cmpi).order-1,Bout,B);
    end
    if nargout > 3 %Not implemented for mvhos yet
        Bpart(:) = interp3(W1,W2,W3,fftshift((me(cmpi).partialbicoh)),P(:),P(:),-P(:)+M(:),'nearest');
        Bpart(isnan(B))=nan;
    end    
    if nargout ==0 %Not implemented for mvhos yet
        ax(cmpi,1)=subplot(1,length(me),cmpi);
        if length(size(Bout))>2
            warning('Plot is averaged over channels')
            while length(size(Bout))>2
                Bout = squeeze(mean(fun(Bout),3));
            end
        end
        imh = pcolor(power,modfreq,fun(Bout));
             set(imh,'facecolor','flat','edgecolor','none');
        if cmpi==1
                ylabel('Envelope modulation freq.(Hz)')
    
         title(ax(cmpi,1),sprintf('Modulogram (Trispectrum Diagonal Slice)'))
        else
         title(ax(cmpi,1),sprintf('Residual Modulogram After Comp. %i',cmpi-1))
        end    
        axis xy
%         xlabel('Envelope modulation freq.(Hz)')
%         ylabel('Band freq. (Hz)')
        xlabel('Band freq. (Hz)')
        
   
    end
end
if nargout >3
    varargout = {Bout,power,modfreq, Bpart };
elseif nargout >0
        varargout = {Bout,power,modfreq};
end