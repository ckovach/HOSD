

function out=filtrecon(x,bsid)

% out = bsident(x,segment,lpfilt,ncomp,varargin)
%
% Reconstruct an input based on the output of bsident. 
%
% See also BSIDENT
%
% C. Kovach 2018

%%% Method to identify impulse
%impulse_method = 'zthreshold';
impulse_method = 'iterz';
thresh = 3;
%impulse_method = 'kmeans';


apply_fun = @(x)x.^3;

xresid = x;

Xrec = zeros(size(x) );
Xfilt = zeros(size(x) );
Ximp = zeros(size(x) );

for k = 1:length(bsid)
%%
      n = size(xresid,1);
     for kk = 1:size(x,2)
        xx = xresid(:,kk);
       

         f = bsid(k).f;
         nX = length(f);
         f(n)=0;
         f = circshift(f,-floor(nX/2));

         h = bsid(k).BFILT;
         h(n) = 0;
         h = circshift(h,-floor(nX/2));

         xfilt = ifft(fft(xx).*fft(h));


        switch impulse_method
            case 'zthreshold'
                ximp = zscore(apply_fun(xfilt))>thresh;
            case 'iterz'
                ximp = isnan(iterz(apply_fun(xfilt),thresh,1));
                
            case 'kmeans'
                kmx = apply_fun(xfilt);
                [km,kmc]  = kmeans(kmx + 0./(zscore(kmx)>0),2);
    %            [km,kmc]  = kmeans(xfilt + 0./(zscore(xfilt)>1),2);
                [~,mxi] = max(kmc);
                ximp = km==mxi;
        end
        xrec = ifft(fft(xfilt.*ximp).*fft(f));%*nX/n;
        a = xrec'*xresid./sum(xrec.^2);
        xrec = a*xrec;
        Xrec(:,kk) = xrec;
        Xfilt(:,kk) = xfilt;
        Ximp(:,kk) = ximp;
    end
     out(k).xrec = Xrec; %#ok<*AGROW>
     out(k).xfilt = Xfilt;
     out(k).ximp = sparse(Ximp);
    out(k).a = a;
    out(k).exvar = 1-sum(abs(xresid-xrec).^2)./sum(abs(xresid).^2);
      xresid =xresid-Xrec;
    
    %%%  A simple measure of compression: 
    %%%     variance explained X total samples / number of values retained 
    %%%     (length of f + number of  impulses used in constructing xrec)
    out(kk).compression = out(kk).exvar*length(x)/(nX+length(out(kk).ximp));
    xresid =xresid-xrec;

end
    
    
