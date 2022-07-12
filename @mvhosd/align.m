 function [Xsh,Xwin,makeplot,iter] = align(me,Xwin,Gpart,maxiter,makeplot,compno)
 % Fit a block of data all at once
% Process input if length is >= buffer size, else add to buffer.
if nargin < 6
    compno = 1;
end
if nargin < 5
    makeplot = true;
end
if nargin < 4 || isempty(maxiter)
    maxiter = 25;
end
%             if ~iscell(xin)
%                 xin = {xin};
%             end
%             nxin = numel(xin);
%             xisnan = isnan(xin);       

if size(Xwin,3)>1
   nsig = size(Xwin,3);
   sigdim = 3;
else
    nsig = size(Xwin,2);
    sigdim = 2;
end

me(1).bufferPos = 0; % Discard the buffer

%             Gpart = me(1).Gpart;

Gpart0 = Gpart;

%      Xwin(end+1:me(1).fftN,:) = 0;
%             Xwin = me(1).Xwin;
Xsh = Xwin;
 Xfilt=nanmean(Xsh,3);


del = Inf;
tol =me(1).sampling_rate/me(1).lowpass(1);
iter = 0;
olddt2 = 0;
olddt =0;

fprintf('\nComponent %3i Iter %3i',compno,0)
color_cycle = 10;
if all(ishandle(makeplot))
    set(makeplot(1:end-1,:),'ydata',zeros(me(1).bufferN,1));
end

if isempty(me(1).sampweight)
    smpw = 1;
else
    smpw = me(1).sampweight;
end

%             std_moment = @(x)nanmean(cumulant(x,me(1).order,1)./(nanmean(x.^2).*nanmean(smpw.^2)).^(me(1).order/2));
%             std_moment = @(x)nanmean(cumulant(x,me(1).order,1)./(nanmean(x.^2)).^(me(1).order/2));
std_moment = @(x)nanmean(cumulant(x,me(1).order,1,true));
switch me(1).order
    case 3
        moment_type = 'skewness';
    case 4
        moment_type = 'kurtosis';
    otherwise
        moment_type = 'standardized cumulant';
end
if me(1).subspace_dim ==0
    nsig = size(Xsh,3);
else
    nsig =  me(1).subspace_dim;
end
if me(1).annealing_start>0
    noise = randn(size(Xsh)).*nanstd(reshape(Xsh,[size(Xsh,1)*size(Xsh,2),1,size(Xsh,3)]))./sqrt(size(Xsh,3));
else
    noise = 0;
end

apa = 1;
while del >tol && iter < maxiter        
    
    if me(1).subspace_dim > 0 && size(me(1).wavefft,3) ~= me(1).subspace_dim
        ff = permute(squeeze(me(1).filterfft)*me(1).projection',[1 3 2]);
    else
        ff = me(1).filterfft;
    end
    [~,plotcompi] = sort(sum(abs(me(1).wavefft).^2.*abs(ff).^2),'descend');

    try
         if any(ishandle(makeplot))
           nplot = size(makeplot,2);
           for k = 1:nplot
                kplot = plotcompi(k);
                if me(1).subspace_dim == 0
                    set(makeplot(1,k),'cdata',Xsh(:,:,kplot)');
                end
                set(makeplot(7,k),'string',sprintf('Ch.%3i, Comp.%3i, Iter.%3i\nMean shift =%2.2fs, %s=%2.2f',kplot,compno,iter,del/me(1).sampling_rate,moment_type,std_moment(Xfilt)));
                set(makeplot(mod(iter,5)+2,k),'ydata',me(1).feature(:,kplot),'Color',hsv2rgb([mod(iter,color_cycle)/color_cycle 1 .8]));
%                             axis tight
                subplot(2,nplot,k+ nplot)
                ylim(minmax(me(1).feature(:))*1.1);
            end
            drawnow
        elseif islogical(makeplot) && makeplot
            figure,
            clear makeplot
            nplot = min(nsig,me(1).maxplotn);
            for k = 1:nplot
                kplot = plotcompi(k);
                subplot(2,max(nplot,me(1).maxplotn),k )
                if me(1).subspace_dim ==0
                    makeplot(1,k) = imagesc(fftshift(me(1).sampt)/me(1).sampling_rate,[],Xsh(:,:,kplot)');
                end
                makeplot(7,k) = title(sprintf('Ch.%3i, Comp.%3i, Iter.%3i\nMean shift =%2.2fs, %s=',kplot,compno,iter,del/me(1).sampling_rate,moment_type ));
                subplot(2,nplot,k+ nplot)
               plh = plot(fftshift(me(1).sampt)./me(1).sampling_rate,me(1).feature(:,kplot)*ones(1,5));
                 for pli = 1:length(plh)

                     set(plh(pli),'Color',hsv2rgb([mod(pli,color_cycle)/color_cycle 1 .8]));
                 end
                 makeplot(2:6,k)=plh;
                 ylim(minmax(me(1).feature(:))*1.1);
            end
%                         xlim([0 me(1).lowpass])

         end
        drawnow
    catch
            plh = false;
    end
    iter=iter+1;
    fprintf('\b\b\b%03i',compno,iter)

    if iter >2
      olddt = me(1).delay;
    end



    if me(1).annealing_start>0
        Xsh0 = Xsh;

%                     noise = randn(size(Xsh))*nanstd(Xsh(:));

        Xsh = Xsh + noise*me(1).annealing_start*me(1).annealing_schedule(iter,maxiter);
        [Xfilt,~,sgn] = me(:,1).apply_mvfilter(Xsh,false,true);
        newdt = me(1).delay;
        FXsh = fft(Xsh0).*conj(fft(me(1).sampt'==newdt)); 
    else
        [Xfilt,FXsh,sgn] = me(:,1).apply_mvfilter(Xsh,false,true);
        newdt = me(1).delay;
    end
    Xsh = real(ifft(FXsh));

%                delt = me(1).radw(me(1).keepfreqs{1})*newdt;
   delt = me(1).radw*newdt;


   %                    
%                if isscalar(smpw)
%                     [~,mxi] = max(Xfilt.^me(1).order);
%                 else
%                     [~,mxi] = max(Xfilt.^me(1).order.*repmat(smpw',size(Xfilt,1),1));
%                 end
%                 if mod(me(1).order,2)==0 
%                     sgn = sign(Xfilt(mxi + (0:size(Xfilt,2)-1)*size(Xfilt,1)).*smpw'); 
%                 end
%                 
%                 samptc=(me(1).sampt);
%                 dt = samptc(mxi);
%                 delt = me(1).radw*dt;
%                 newdt = dt;
   %  deltg = delt(me(1).keepfreqs{1},:);
    Gpart = (sgn.*exp(-1i.*delt)).*Gpart; 
    G = nanmean(Gpart,2);

%                Xsh = real(ifftshift(ifft(exp(1i*delt).*FXwin),1));
%                Xsh = real(ifft(exp(1i*delt).*FXwin));
    features = nanmean(Xsh,2);

    %%% G should be matched to features so as to produce a
    %%% peak at zero lag. Might as well enforce this
    %%% explicitly as aligment arrors might accrue.

    L = ifft(G.*fft(features));
    [~,mxi] = max(real(L).^(2-mod(me(1).order,2)));

    if mod(me(1).order,2)==0
        sgc=sign(real(L(mxi+permute((0:size(L,3)-1)*size(L,1),[1 3 2]))));
    else
        sgc=1;
    end
    flcorrection = me(1).sampt(mxi);
    delt2 = me(1).radw*flcorrection; 
    Gpart = Gpart.*exp(1i*permute(delt2,[1 3 2])).*sgc;

    if size(G,1) == me(1).fftN
        G = G(me(1).keepfreqs{1},:,:);
    end

     if me(1).subspace_dim > 0
%                     [u,l,v] = svds(squeeze(features),me(1).subspace_dim);
        [u,l,v] = svds(squeeze(G),me(1).subspace_dim);
%         me(1).feature = permute(squeeze(features)*real(v),[1 3 2]);
        me(1).waveftlag =fft(features); 
        me(1).G = permute(squeeze(G)*real(v),[1 3 2]);
        me(1).projection = real(v);
    else

        me(1).G = G;%.*apa;
        me(1).waveftlag= fft(features);

     end           


    if me(1).adjust_lag
       if me(1).subspace_dim > 0 && size(me(1).filterftlag,3) == size(me(1).projection,2)
           ff = permute(squeeze(me(1).filterftlag)*me(1).projection',[1 3 2]);
       else
           ff = me(1).filterftlag;
       end
       ffun = ifftshift(sum(real(ifft(ff.*abs(me(1).waveftlag+eps))),3),1);                   
       mph = sum(exp(-1i*2*pi*me(1).sampt(:)./me(1).fftN).*sum(abs(ffun).^2,2))./sum(sum(sum(abs(ffun).^2,2),3));                   
       mph = mph./(abs(mph)+eps);
       me(1).lag = mph; % Circularshift to keep filter energy centered on the window
    end
%                 Gpart = Gpart;
    Xwin = Xsh;


    % checks two and one step back to reduce getting
    % trapped at points of cyclical stability.
    del = std(olddt2-newdt);%min(sqrt(mean((olddt-newdt).^2)),sqrt(mean((olddt2-newdt).^2)));
    olddt2 = olddt;



end

%%% Set the delays to the correct value for the original
%%% data set;
 [~,~] = me(1).apply_mvfilter(Xwin,false,true);
%             T = T+ repmat(me(1).delay,size(T,1),1);
%             T(T<1)=1;
%             T(T>length(xin))=length(xin);
if size(me,2)>1
   Xrec = me(:,1).xrec(Xwin);
%                if nargout > 0
%                    [Xsh2,Xwin2,T2] = me(2:end).get_block(xin-xrec,maxiter,makeplot,compno+1);
%                    me(2).Xwin = Xwin-Xrec;
       delt = me(1).radw*me(1).delay; 
%                    me(2).Gpart = (Gpart - me(1).filterfft).*exp(1i.*delt); % The partial delay filter has to retain the correct alignment to Xwin.
%                    [Xsh2,Xwin2] = me(2:end).get_block(xin-xrec,maxiter,makeplot,compno+1);
       [Xsh2,Xwin2] = me(2:end).align(Xwin-Xrec,(Gpart - me(1).filterfft).*exp(1i.*delt),maxiter,makeplot,compno+1);
       Xsh = cat(3,Xsh,Xsh2);
       Xwin = cat(3,Xwin,Xwin2);
%                    T = cat(3,T,T2);
%                else
%                     me(2:end).get_block(xin-xrec,maxiter,makeplot,compno+1);
%                end
end
end
