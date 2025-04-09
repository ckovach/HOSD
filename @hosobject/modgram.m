function varargout = modgram(me,Bplot0)

%Extract the modulogram from a 4th order spectrum

if me(1).order~=4
    error('Spectrum is not 4th order');
end

fun = @(x)x;

for cmpi = 1:length(me)
    if nargin < 2
        Bplot = fftshift(me(cmpi).bicoh);
        fun = @abs;
    elseif isa(Bplot0,'function_handle')
        Bplot = fftshift(me(cmpi).bicoh);
        fun = Bplot0;
    elseif size(Bplot0,2) == length(me)
        Bplot = Bplot0(:,cmpi);
        Bplot = fftshift(Bplot(me(cmpi).fullmap));
    else 
        Bplot = Bplot0;
    end

    wb0 = cellfun(@fftshift,me(cmpi).freqindx.Bfreqs,'uniformoutput',false);
         
    power = wb0{1};
   % power = power(power>0 & power >me(cmpi).highpass &power< me(cmpi).lowpass);
    power = power(power>0);
    modfreq = wb0{end}-min(wb0{end});
  %  modfreq = modfreq(abs(modfreq)>me(cmpi).shighpass & abs(modfreq)<=me(cmpi).slowpass);
    modfreq = modfreq(abs(modfreq)>0);
    
    [P,M] = meshgrid(power,modfreq);
    [W1,W2,W3] = meshgrid(wb0{:});
    B = nan*M;
    Bpart = B;
    B(:) = interp3(W1,W2,W3,Bplot,P(:),P(:),-P(:)+M(:));
    if nargout > 3
        Bpart(:) = interp3(W1,W2,W3,fftshift((me(cmpi).partialbicoh)),P(:),P(:),-P(:)+M(:));
        Bpart(isnan(B))=nan;
    end    
    if nargout ==0
        ax(cmpi,1)=subplot(1,length(me),cmpi);
        imh = pcolor(power,modfreq,fun(B));
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
    varargout = {B,power,modfreq, Bpart };
elseif nargout >0
        varargout = {B,power,modfreq};
end