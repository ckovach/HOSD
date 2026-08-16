function varargout=reseed(s)

%
% RESEED seeds the uniform and normal random generators, RAND and RANDN
% with the current time obtained with NOW.
%
% See also: RAND, RANDN, NOW.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

%For legacy reseeding
 umethod = 'twister';
 nmethod = 'state';
% 
if nargin < 1 || isempty(s)
    s = typecast(now,'uint32');
    s=s(1); %Using the least significant 4 bytes of the value returned by now
end

try
    rng(s)
catch
    rand(umethod, double( s )); 
    srandn=double( swapbytes( s ) );
    randn( nmethod, srandn ); %using swapbytes so that the randn and rand seeds are
                                             %  at least somewhat arbitrary
                                             %with respect to each other.
end

if nargout > 0
    varargout{1}=s;
end

