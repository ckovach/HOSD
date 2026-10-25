function varargout =plotpolycoh(hos,ax,plotwhat,absfun)

%Plot polycoherence. For orders greater than 3, plot 2d slices.
%
% plotpolycoh(hos,ax,plotwhat,absfun)
%   hos      - a hosobject, or an array of them (one plot each)
%   ax       - the axes to draw into (default: subplots of the current figure)
%   plotwhat - what is drawn: 'abs' (the default; also 'mag', 'magnitude'),
%              'phase' (also 'angle'), a numeric scalar t (the phase where
%              the magnitude exceeds t, NaN elsewhere; 'phase' is t = 0) or
%              a function handle.
%   absfun   - a function handle applied instead of the one plotwhat names
%              (older callers pass the same handle as plotwhat and absfun);
%              empty means none, and any other value is refused.
%
% The values drawn are looked up directly in the reduced coefficient vector
% through the square-form map (hosobject/square_index, freqindx.remap and
% PDconj) instead of being read out of the full square-form array through
% interp2 / interp3 / interpn: the modulogram lattice (p,p,m-p), the
% bispectral diagonal (f,f) and the 2-D slices of a 4th-order spectrum all
% lie on the frequency grid, so the nearest-neighbor interpolation was an
% exact lookup at the cost of forming, shifting and interpolating arrays of
% size prod(dims). The plotting function is applied once to the looked-up
% values in every branch, the NaN of points off the grid included, so it is
% expected to act elementwise (abs, real, angle, a threshold ...).

if nargin < 2 
    ax = [];
end
% Axes may be passed as graphics objects or as numeric handles (older
% callers); numeric handles are converted so that property access and
% indexed assignment below work the same way for both.
if ~isempty(ax) && isnumeric(ax)
    ax = handle(ax);
end
% Every graphics call below targets the axes it draws into. When the
% caller supplies the axes (e.g. a uiaxes inside a uifigure) nothing may
% fall through to gca / gcf: with no classic figure open those calls
% would create a new empty figure titled 'Component 1'.
ownAxes = isempty(ax);
thresh = 0;
if nargin < 3 || isempty(plotwhat)
    plotwhat = 'abs';
elseif isnumeric(plotwhat) && isscalar(plotwhat)
    thresh = plotwhat;
    plotwhat = 'phase';
end
if isa(plotwhat,'function_handle')
    fun = plotwhat;
else
    name = '';      % anything that is not text reaches the error below
    if (ischar(plotwhat) && isrow(plotwhat)) || (isstring(plotwhat) && isscalar(plotwhat) && ~ismissing(plotwhat))
        name = lower(char(plotwhat));
    end
    switch name
        case {'mag','abs','magnitude'}
            fun = @abs;
        case {'phase','angle'}
            fun = @(x)angle(x) + 0./(abs(x)>thresh);
        otherwise
            error('plotpolycoh:plotwhat', ...
                'Unrecognized plotting option: plotwhat must be ''abs'', ''phase'', a numeric threshold or a function handle.')
    end
end
if nargin < 4 || isempty(absfun)
    absfun = fun;
elseif ~isa(absfun,'function_handle')
    error('plotpolycoh:absfun','absfun (the 4th argument) must be a function handle, or empty.')
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
subxy = almostSquare(length(hos));
txth = [];
imh = gobjects(0);      % an object array, so imh(k) = imagesc(...) keeps graphics handles
for hi = 1:length(hos)
   
    h = hos(hi);
    fx = h.freqindx;
    if h.diagonal_slice && h.order == 4
        subxy = [length(hos) 2];
        wb0 = cellfun(@fftshift,fx.Bfreqs,'uniformoutput',false);
     
        power = wb0{1};
        power = power(power>0 & power >h.highpass &power< h.lowpass);
        modfreq = wb0{end}-min(wb0{end});
        modfreq = modfreq(abs(modfreq)>h.shighpass & abs(modfreq)<=h.slowpass);
        
        [P,M] = meshgrid(power,modfreq);
        B = lookup(h,[P(:),P(:),-P(:)+M(:)],h.bicohreduced,true);
        Bpart = lookup(h,[P(:),P(:),-P(:)+M(:)],h.partialbicohreduced,true);
        B = reshape(B,size(P));
        Bpart = reshape(Bpart,size(P));
        Bpart(isnan(B))=nan;
         if hi>size(ax,1)
             ax(hi,1) = subplot(subxy(1),subxy(2),1 + 2*(hi-1));
             if size(ax,2)>1
                ax(hi,2) = subplot(subxy(1),subxy(2),2 + 2*(hi-1));
             end
%         else
%             axes(ax(hi))
         end
        % double: pcolor refuses logical CData (a handle such as @(x)abs(x)>t)
        imh = pcolor(power,modfreq,double(absfun(B)),'parent',ax(hi,1));
        set(imh,'facecolor','flat','edgecolor','none');
        if hi==1
         title(ax(hi,1),sprintf('Modulogram (Trispectrum Diagonal Slice)'))
        else
         title(ax(hi,1),sprintf('Residual Modulogram After Comp. %i',hi-1))
        end    
        axis(ax(hi,1),'xy')
%         xlabel('Envelope modulation freq.(Hz)')
%         ylabel('Band freq. (Hz)')
        ylabel(ax(hi,1),'Envelope modulation freq.(Hz)')
        xlabel(ax(hi,1),'Band freq. (Hz)')
        
         if hi==1
            cax = caxis(ax(hi,1));
        else
            caxis(ax(hi,1),cax)
     
         end
         if size(ax,2)>1
            imh = pcolor(power,modfreq,double(absfun(Bpart)),'parent',ax(hi,2));
            set(imh,'facecolor','flat','edgecolor','none');
            title(ax(hi,2),sprintf('Component %i Trispectrum Modulogram',hi))
            axis(ax(hi,2),'xy')
            caxis(ax(hi,2),cax)
         end
    elseif h.diagonal_slice && h.order == 3
        wb0 = cellfun(@fftshift,fx.Bfreqs,'uniformoutput',false);
        f = wb0{1};
        B = nan*f;
        Bpart = B;
        B(:) = absfun(lookup(h,[f(:),f(:)],h.bicohreduced,true));
        Bpart(:) = absfun(lookup(h,[f(:),f(:)],h.partialbicohreduced,true));
         if hi>length(ax)
             ax(hi) = subplot(subxy(1),subxy(2),hi);
%         else
%             axes(ax(hi))
         end
        imh = plot(f,B.*(f>=0)+ (1-(f>=0)).*Bpart,'parent',ax(hi));
         title(ax(hi),sprintf('Component %i Bispectrum Diagonal Slice',hi))
         xlabel(ax(hi),'freq.(Hz)')
         if hi==1
             cax = caxis(ax(hi));
         else
             caxis(ax(hi),cax);
         end
    else
        wb0 = cellfun(@fftshift,fx.Bfreqs,'uniformoutput',false);
        switch scale

            case 'log'
                wb = cellfun(@(x)log10(logspace(log10(h.highpass),log10(max(x)),length(x))),wb0,'uniformoutput',false);
                do_interp = true;
  
            case {'lin','linear'}
                wb = cellfun(@(x)x(x>h.highpass),wb0,'uniformoutput',false);
                do_interp = true;
                
            otherwise
                wb=wb0;
                do_interp=false;
  
        end
      
        wb(end+1) = {0};

%         ax(hi) = subplot(subxy(1),subxy(2),hi,'UserData',wb);
        if hi>length(ax)
             ax(hi) = subplot(subxy(1),subxy(2),hi,'UserData',wb);
%         else
%             axes(ax(hi))
         end
        hold(ax(hi),'on')

        dwb = cellfun(@(x)diff(x([1 end])),wb(1:2));

        nsamps = cellfun(@length,wb);

        %%% Slices along the axes beyond the first two: one image per
        %%% combination (half of the last axis, the rest being redundant
        %%% under conjugation).
        if length(wb)>3
            K = cellfun(@(x)(1:length(x))',wb(3:end-1),'uniformoutput',false);
            if isscalar(K)
                K = K{1};
            else
                [K{:}] = ndgrid(K{:});
                K = cellfun(@(x)x(:),K,'uniformoutput',false);
                K = [K{:}];
            end
            K = K(K(:,end)<=length(wb{end-1})/2,:);
        else
            K = 1;
        end
        dims = [];
        for k = 3:length(wb)
            dims(k-2,:) = almostSquare(nsamps(k)); %#ok<AGROW>
        end    

        [W1,W2] = ndgrid(wb{1},wb{2});
        if do_interp
            %%% Resampled axes: interpolate the full array (rare; the scale
            %%% is fixed to the grid above).
            wbin = wb0; [wbin{:}] = ndgrid(wb0{:});
            wbout = wb(1:end-1); [wbout{:}] = ndgrid(wb{1:end-1});
            Bfull = interpn(wbin{:},absfun(fftshift(h.bicoh)),wbout{:});
        end
        for k = 1:size(K,1)

            cent = [0 0];
            for kk = 1:size(K,2)
                cent =  cent + [mod(K(k,kk)-1,dims(kk,2)), floor((K(k,kk)-1)/dims(kk,2))].*prod(dims(1:kk-1,:),1);
            end        
            wcent = cent.*dwb;

            imk = num2cell(K(k,:));
            if do_interp
                Bslice = Bfull(:,:,imk{:});
            else
                %%% The slice's frequencies on the axes beyond the second
                wrest = zeros(1,length(wb)-3);
                for kk = 1:length(wrest)
                    wrest(kk) = wb{kk+2}(K(k,kk));
                end
                Wq = [W1(:),W2(:),repmat(wrest,numel(W1),1)];
                Bslice = reshape(absfun(lookup(h,Wq,h.bicohreduced,true)),size(W1));
            end
            if plot_partialbc && ~do_interp
                Wq = [W1(:),W2(:),repmat(wrest,numel(W1),1)];
                Bpslice = reshape(absfun(lookup(h,Wq,h.partialbicohreduced,true)),size(W1));
            else
                Bpslice = Bslice;
            end
            Bplot = Bslice.*(W2>=W1) + (1-(W2>=W1)).*Bpslice;
            imh(imk{:}) = imagesc(wb{1}+wcent(1),wb{2}+wcent(2),Bplot,'parent',ax(hi));
            hIm = imh(imk{:});
            if isnumeric(hIm)                     % a numeric handle array: convert for property access
                hIm = handle(hIm);
            end
            delete(datatip(hIm));
%             delete(dt)
            Ws = {W1,W2};
            for kk = 3:length(wb)-1
                Ws{kk} = wb{kk}(K(k,kk-2))*ones(size(W1));
            end
            Ws{length(wb)} = -sum(cat(3,Ws{:}),3);
            for kk = 1:length(Ws)
                q=hIm.DataTipTemplate;
                q.DataTipRows(kk) = dataTipTextRow(sprintf('f_%i:',kk),Ws{kk});
            end
             q=hIm.DataTipTemplate;
             q.DataTipRows(length(Ws)+1) = dataTipTextRow('Value:',double(Bslice));
            txt = {};
            for kk = 3:length(wb)-1
                if kk==3
                    txt{kk-2} = sprintf('%0.2f',wb{kk}(K(k,kk-2))); %#ok<AGROW>
                else
                    txt{kk-2} = sprintf(', %0.2f',wb{kk}(K(k,kk-2))); %#ok<AGROW>
                end
            end
            if ~isempty(txt)
             txth(imk{:})=text(wcent(1)+dwb(1)/2+wb{1}(1),wcent(2)+.9*dwb(2)+wb{2}(1),[txt{:}],'Color','w','HorizontalAlignment','center','fontsize',8,'parent',ax(hi)); %#ok<AGROW>
            end
        end
        title(ax(hi),sprintf('Component %i',hi))
        axis(ax(hi),'image')
        if hi==1
            cax = caxis(ax(hi));
        else
            caxis(ax(hi),cax)
        end
        if nargout >= 3
            %%% The third output of this branch has always been the whole
            %%% (shifted) polycoherence array; it is formed only on request.
            if do_interp
                B = Bfull;
            else
                B = absfun(fftshift(h.bicoh));
            end
        end
%         tick = @(x)round(linspace(x(1),x(end),10));
%          set(gca,'xtick',tick(wb{1}),'ytick',tick(wb{2}),'xticklabel',ticklbl(wb{1}),'yticklabel',ticklbl(wb{2}))
    end
    if nargout >=3
        Bs{hi} = B; %#ok<AGROW>
    end
end

if ownAxes
    colorbar(ax(end)) % a caller-supplied axes manages its own colour scale
end

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
% set(gca,'xticklabel',lblx(get(ax,'xtick')),'yticklabel',lblx(get(ax,'ytick')));
%     
    
function vals = lookup(h,W,v,conjugate)

% Coefficients of the reduced-form vector v at the square-form points whose
% first order-1 frequencies are the rows of W: v(remap), conjugated where
% PDconj is set (as hos.bicoh does) when conjugate is true, NaN off the
% grid.

fx = h.freqindx;
[lin,inside] = h.square_index(W);
vals = nan(size(W,1),1);
j = double(fx.remap(lin(inside)));
b = v(j);
if conjugate
    cj = fx.PDconj(lin(inside));
    b(cj) = conj(b(cj));
end
vals(inside) = b;

function out =  almostSquare(N)

rg = 0:sqrt(N);
divs = arrayfun(@(x)find(rem(x,1:x)==0),N+rg,'uniformoutput',false);
divr = arrayfun(@(x,y)y{1}-x./y{1},N+rg,divs,'uniformoutput',false);

[mn,mni] = cellfun(@(x,y)min(abs(x)),divr,divs);
[~,mni2] = min(mn);
dm = divs{mni2}(mni(mni2));
out = [(N+rg(mni2))/dm  dm];
