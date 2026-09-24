
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
% Iout(i,j) is the number of elements of mapto with subscript i along dim
% and value j. It is built directly from the nonzero, unmasked elements
% with one call to sparse (which accumulates repeated (i,j) pairs): one
% pass over the array to find the nonzeros, then O(nnz log nnz) work on
% those, never a matrix the size of numel(mapto). The earlier form
% composed a numel(mapto) x veclength map with a numel(mapto) x
% size(mapto,dim) kron chain, (Imap'*Isum)', an order of magnitude slower
% on a 307^3 square form.
%
% See also FREQ2INDEX

%
% C. Kovach 2019, 2026
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

%%% Elements that contribute: nonzero (one pass over the array) and not
%%% masked (checked on the nonzero values only).
fmapto = find(mapto);
vals = double(mapto(fmapto));
maskval = maskval(maskval ~= 0);
if ~isempty(maskval)
    drop = false(size(vals));
    for v = maskval(:)'
        drop = drop | vals == v;
    end
    fmapto = fmapto(~drop);
    vals = vals(~drop);
end

%%% Subscript along dim of each contributing element, without ind2sub
%%% (only the one subscript is needed).
stride = prod(sz(1:dim-1));
subdim = mod(floor((fmapto-1)/stride),sz(dim))+1;

%%% Iout(i,j) counts the elements with subscript i along dim that map to j.
Iout = sparse(subdim,vals,1,sz(dim),veclength);
