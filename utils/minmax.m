% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

function [mnmx,mnmxi] = minmax(X)

if nargout <2
    mnmx = cat(1,min(X),max(X));
else
   [mn,mni] = min(X);
   [mx,mxi] = max(X);
   mnmx = [mn;mx];
   mnmxi = [mni;mxi];
end