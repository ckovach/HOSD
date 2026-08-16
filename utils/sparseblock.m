
function SpX = sparseblock(X,b, varargin)

%
% spbX = sparseblock(X,b)
%
%   SPARSEBLOCK creates a block-diagonal sparse matrix from X divided 
%   according to b: if X is an M x N matrix, then each block, Bi, has 
%   dimension M x b(i), or M x b if b is a scalar which evenly divides N. 
%   For example if
%
%                 >> X =   [ 1 2 3 ; 4 5 6 ]
%
%              
%                   X =
%   
%                        1     2     3
%                        4     5     6
%   then
%         
%                 >> Y = full( sparseblock(X,[2 1]) )
%                 
%                 Y =
% 
%                          1     2     0
%                          4     5     0
%                          0     0     3
%                          0     0     6
%
% 
%   
% spbX = sparseblock(X,b,'transpose') 
% 
%    Same as above, except that X and SpX are transposed so that the expansion 
%    occurs across colums rather than rows.
%
%   
%See also UNSPARSIFY

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

%
% Written by C. Kovach 2007
%

transpose = false;

i=1;
while i <= length(varargin)
    switch lower(varargin{i})
        case 'transpose' %Transposes the input and the output if true
            transpose = true;       
        otherwise
            error('%s is not a valid keyword.',varargin{i})
    end
    i = i+1;
end

if ~isa(X,'double')
    error('Input must be class double.');
end

if isempty(X) || isempty(b)
    SpX = sparse(X);
    return
end
    
if transpose
    X = X';
end

if length(size(b)) > 2 || ~any(size(b) == 1)
    error('Second argument must be a vector.');
end
    
b = b(:);
nblocks = length(b);

if nblocks == 1
    nblocks = size(X,2)./b;
    if mod(nblocks,1) ~= 0
        error('b doesn''t evenly divide the number of Columns of X.')
    end
    b = b*ones(nblocks,1);
end

Jc = (0:size(X,1):numel(X));


strows = zeros(size(X,2),1);
strows( cumsum(b(1:end-1))+1) = 1;
strows = cumsum(strows)*size(X,1);

Ir = kron(strows,ones(size(X,1),1)) + repmat((1:size(X,1))', size(X,2),1)-1;

nrows = size(X,1)*nblocks;
SpX = sparseblockmex(X,Ir,Jc, nrows ); %sparseblockmex is a mex function

if transpose
    SpX = SpX';
end

