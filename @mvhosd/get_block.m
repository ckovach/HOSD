 function [A,B,makeplot] = get_block(me,in,maxiter,makeplot,segment,compno,initialize)

      %Process a block of data. 
      
      if nargin < 6 || isempty(compno)
          compno = 1;
      end
      if nargin < 5 || isempty(segment)
          segment = [];
      end
      if nargin < 7 || isempty(initialize)
          initialize = true;
      end
      if nargin < 4 || isempty(makeplot)
          makeplot = true;
      end
      if nargin < 3 || isempty(maxiter)
          maxiter = 25;
      end
      if size(in,2) > 1 && size(in,3) ==1
           in = permute(in,[1 3 2]);
      end

      if all(isnan(in))
          A = [];
          B = [];
          return
      end

        nfp = 0;
        nchan = size(in,3);
        for k = 1:nchan
            nfp=fprintf([repmat('\b',1,nfp),'\nComp. %i,estimating HOS for chan. %i'],compno,k)-nfp;
            X(:,:,k) = me(1).chop_input(in(:,k));
            Gpart(:,:,k) = me(1).partial_delay_filt(X(:,:,k),true,true); 
            if k==1
                X(:,:,nchan)=0;
                Gpart(:,:,nchan)=0;
            end
        end


        me(1).use_adaptive_threshold=false;

        feature = nanmean(X,2);
        filterfft = nanmean(Gpart,2); %#ok<*NANMEAN>
        if me(1).subspace_dim > 0
%                     [u,l,v] = svds(squeeze(feature),me(1).subspace_dim);
            [u,l,v] = svds(squeeze(filterfft),me(1).subspace_dim);
%             me(1).feature = permute(squeeze(feature)*real(v),[1 3 2]);
            me(1).feature = permute(squeeze(feature),[1 3 2]);
            me(1).filterfft = permute(squeeze(filterfft)*v,[1 3 2]);
            me(1).projection = real(v);
        else
            me(1).feature = feature;
            me(1).filterfft = filterfft;
            me(1).projection = 1;
        end
%                 if me(1).do_static_ica >0
%                     for rep = 1:10
%                        [Xfilt,FXsh,sgn,mvXfilt] = me(:,1).apply_mvfilter(X,false,true);
%                        D = real(ifft(exp(1i*me(1).radw*me(1).delay)).*sgn);
%                        XFsh = real(ifft(fft(mvXfilt).*fft(D)));
%                      %  XFsh = XFsh(1,:,:);
%                        XFsh = XFsh(abs(me(1).sampt)<=ceil(me(1).buffersize/4),:,:);
%     %                    XFsh = XFsh./sqrt(nansum(nansum(XFsh.^2,3),2)); %Normalizing total energy in each window in order to reduce the tendency for the spatial weighting to collapse onto outliers
%                        mvxfrs = reshape(XFsh,size(XFsh,1)*size(XFsh,2),size(XFsh,3));
%                        Apica = pica(mvxfrs(~any(isnan(mvxfrs),2),:),1,me(1).order,ones(size(mvxfrs,2),1));
%                        me(1).G = me(1).G.*permute(Apica,[2 3 1]);
%                     end
%                 end

        if ~me(1).do_static_ica
            [A,B,makeplot] = me(1).align(X,Gpart,maxiter,makeplot,compno);
        else
          niter = 0;
          while niter <= maxiter
             [A,B,makeplot,iter] = me(1).align(X,Gpart,min(maxiter,me(1).do_static_ica),makeplot,compno);
             niter = niter+me(1).do_static_ica;


             apvar=Inf;
             rep = 0;
%                  apa = 1;
            AA = 1;
             while apvar>1e-3 && rep <= 10
                   [xf,XF] = me(:,1).xfilt(A);
                   XF = XF(abs(me(1).sampt)<=me(1).buffersize/4,:,:);
                   XF = XF(:,~any(isnan(sum(XF,3)),1),:);
                   XF = reshape(XF,size(XF,1)*size(XF,2),size(XF,3));
                   Apica = pica(XF,1,me(1).order,ones(size(XF,2),1),[],false);
                   if me(1).subspace_dim ==0
                       me(1).G = me(1).G.*permute(Apica,[2 3 1]);
                   else
                       me(1).projection = Apica.*me(1).projection;
                   end
                   AA = AA.*Apica;

%                        apvar =var(Apica); 
                   apvar = mean((Apica-1).^2);%This converges to 1
%                        apa = apa.*permute(Apica,[2 3 1]);
                   rep=rep+1;
             end
%                      if me(1).subspace_dim ==0
%                          Gpart = Gpart.*permute(AA,[2 3 1]);
%                      end
             [~,FXsh] = me(:,1).apply_mvfilter(X,false,true);

              Xsh = real(ifft(FXsh));
              if me(1).subspace_dim == 0
                  me(1).feature = nanmean(Xsh,2);
                  Gpart = Gpart.*permute(AA,[2 3 1]);
              else
                  feature = squeeze(nanmean(Xsh,2));
                  delt = me(1).radw*me(1).delay;     
                  filterfft = nanmean((exp(-1i.*delt)).*Gpart,2);
%                           [v,l] = svds(me(1).projection,me.subspace_dim);
%                          [u,l,v] = svds(feature,me(1).subspace_dim);
%                              wgtfun = ifft(fft(feature).*abs(squeeze(filterfft)));
%                             [u,l,v] = svds(wgtfun,me(1).subspace_dim);                 
                   [u,l,v] = svds(squeeze(filterfft),me(1).subspace_dim);
                  me(1).feature = permute(feature,[1 3 2]);
%                   me(1).feature = permute(feature*real(v),[1 3 2]); 
                   me(1).projection = real(v);

%                    me(1).filterfft = permute(squeeze(filterfft)*v,[1 3 2]);
                   me(1).filterfft = permute(squeeze(filterfft)*v,[1 3 2]);
              end
              A = Xsh;
              if rep <2  && iter <3
                     break
              end
           end
        end
        if size(A,4) ==1
            A = permute(A,[1 2 4 3]);
        end
%                  [~,xfilt] = me(1).xfilt(permute(zresid,[1 3 2]));

%                 if nk > 1 && skf(nk-1)< params.skewness_threshold && skf(nk)< params.skewness_threshold 
%                     fprintf('Skewness under %0.3f for the last 2 components... stopping at %i.',params.skewness_threshold,nk)
%                     return
%                 else
        if length(me)>1
            xrec = me(1).xrec(in);
             [A2,B2] = me(2:end).get_block(in-xrec,maxiter,makeplot,segment,compno+1,initialize);
            A = cat(3,A,A2);
            B = cat(3,B,B2);
        end

  end