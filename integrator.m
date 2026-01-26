
function Iout = integrator(mapto,dim,veclength,maskval)

%
% Iout = integrator(mapto,dim,veclength,maskval)
%
% Create a sparse matrix that integrates over all but 1 dimensions of a K-dimensional
% array whose elements map to those of a vector according to mapto. Integration is over
% all dimensions except dim.
%
% Outputs:
%     Iout: A matrix of dimension size(mapto,dim) x veclength matrix. 
%           Given a column vector, c, for which D = c(mapto):
%           Iout*c = sum(...sum(sum( D , d_1 ),d_2), ... d_{K-1})
%               where d_k ~= dim.
%           This avoids the need for an explicit reconstruction of D given
%           the non-redundant representation of the elements in D as c.
% Inputs:
%     mapto: an array of the same size as the target array, whose
%                 elements are indices into a vector. 
%     dim:   Dimension to preserve (default = 1).
%     veclength: Length of the source vector (default == max(mapto(:))).
%     maskval: Values in mapto corresponding to entries in D to ignore, 
%              assumed to be zero-valued (defualt = 0).
%
% See also FREQ2INDEX

%
% C. Kovach 2019
%

sz = size(mapto);

if nargin < 2 || isempty(dim)
    dim = 1;  
end
if nargin < 3 || isempty(veclength)
    veclength = double(max(mapto(:)));
end
if nargin < 4 
    maskval = 0;
end

if ~isequal(maskval,0)
    mapto(ismember(mapto,maskval))=0;
end

inds = arrayfun(@(x)1:x,sz,'uniformoutput',false);

% Is = repmat({[]},1,length(sz));

% [Is{:}] = ndgrid(inds{:});

%%% This matrix maps elements in the vector to elements in the array; 
%%% that is, reshape(x'*Imap,size(mapto)) casts x into the array according 
%%% to mapto.
fmapto = find(mapto);
Imap = sparse(fmapto,double(mapto(fmapto)),ones(size(fmapto)),numel(mapto),veclength);


%%% This matrix integrates over all dimensions in the array but dim
Isum = sparse(1);
for k = 1:length(inds)
    if k == dim
        Isum = kron(speye(sz(k)),Isum);
    else
        Isum = kron(ones(sz(k),1),Isum);
    end
end
%     %Very slow way of doing it
%     Isum = arrayfun(@(x)sparse(Is{dim}(:)==x),inds{dim},'uniformoutput',false);
%     Isum = [Isum{:}]; 

%%% This matrix combines the two operations...
Iout = (Imap'*Isum)';


