function varargout = modgram(me,Bplot0)

%Extract the modulogram from a 4th order spectrum
%
%  [Bout,power,modfreq,Bpart] = modgram(hos)       polycoherence; plotted as |.|
%  [...] = modgram(hos,fun)                        fun applied to the plotted values
%  [...] = modgram(hos,v)                          a reduced-form vector the length of
%                                                  hos.B (one column per component
%                                                  for an array), e.g. a mask
%
% The modulogram is the diagonal slice of the trispectrum read along the
% modulation axis: its value at carrier p and modulation frequency m is the
% coefficient at (w1,w2,w3) = (p,p,m-p), w4 = -(p+m). Those points lie on
% the frequency grid, so they are looked up directly in the reduced
% coefficient vector through freqindx.remap (see SQUARE_INDEX), conjugated
% where freqindx.PDconj says so as hos.bicoh does, and NaN where the point
% is off the grid or was not estimated -- the values the earlier
% nearest-neighbor interp3 over the full square-form array returned,
% without forming that array. A supplied vector is read through the map as
% is (no conjugation), as before. Bout is modulation x carrier (x channel).

if me(1).order~=4
    error('Spectrum is not 4th order');
end

fun = @(x)x;
usevector = false;
if nargin < 2
    fun = @abs;
elseif isa(Bplot0,'function_handle')
    fun = Bplot0;
else
    usevector = true;
end

for cmpi = 1:length(me)
    h = me(cmpi);
    fx = h.freqindx;
    len = size(fx.Is,1);

    wb0 = cellfun(@fftshift,fx.Bfreqs,'uniformoutput',false);
    power = wb0{1};
    power = power(power>0);
    modfreq = wb0{end}-min(wb0{end});
    modfreq = modfreq(abs(modfreq)<h.slowpass);
    modfreq = modfreq(abs(modfreq)>0);
    [P,M] = meshgrid(power,modfreq);

    [lin,inside] = h.square_index([P(:),P(:),-P(:)+M(:)]);
    j = ones(size(lin));
    j(inside) = double(fx.remap(lin(inside)));
    cj = false(size(lin));
    cj(inside) = fx.PDconj(lin(inside));

    if usevector
        if size(Bplot0,2) >= cmpi && size(Bplot0,2) == length(me)
            v = Bplot0(:,cmpi);
        else
            v = Bplot0(:,1);
        end
        if numel(v) ~= len+1
            error('The vector must have length(hos.B) = %d entries, one per reduced coefficient plus the trailing slot.',len+1);
        end
        vals = {v};                                  % read as is, as Bplot0(fullmap) did
        conjugate = false;
    else
        BC = h.bicohreduced;                         % (len+1) x 1 x channels
        vals = cell(1,size(BC,3));
        for k = 1:size(BC,3)
            vals{k} = BC(:,1,k);
        end
        conjugate = true;
    end

    Bout = [];
    for k = 1:length(vals)
        B = nan(size(P),'like',double(vals{k}(1)));
        v = vals{k};
        b = v(j(inside));
        if conjugate
            b(cj(inside)) = conj(b(cj(inside)));
        end
        B(inside) = b;
        Bout = cat(3,Bout,B);
    end

    Bpart = nan(size(P));
    if nargout > 3 %Not implemented for mvhos yet
        pBC = h.partialbicohreduced;
        b = pBC(j(inside));
        b(cj(inside)) = conj(b(cj(inside)));
        Bpart(inside) = b;
        Bpart(isnan(B))=nan;
    end
    if nargout ==0 %Not implemented for mvhos yet
        ax(cmpi,1)=subplot(1,length(me),cmpi); %#ok<AGROW>
        Bplot = Bout;
        if length(size(Bplot))>2
            warning('Plot is averaged over channels')
            while length(size(Bplot))>2
                Bplot = squeeze(mean(fun(Bplot),3));
            end
        else
            Bplot = fun(Bplot);
        end
        imh = pcolor(power,modfreq,Bplot);
             set(imh,'facecolor','flat','edgecolor','none');
        if cmpi==1
                ylabel('Envelope modulation freq.(Hz)')
    
         title(ax(cmpi,1),sprintf('Modulogram (Trispectrum Diagonal Slice)'))
        else
         title(ax(cmpi,1),sprintf('Residual Modulogram After Comp. %i',cmpi-1))
        end    
        axis xy
        xlabel('Band freq. (Hz)')
    end
end
if nargout >3
    varargout = {Bout,power,modfreq, Bpart };
elseif nargout >0
        varargout = {Bout,power,modfreq};
end
