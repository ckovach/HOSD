
function [X,b] = unsparsify( SpX , varargin)

%
% [X,b] = unsparsify( SpX )
%
%     From a sparse block diagonal matrix, SpX, where the ith block, Bi, 
%     has size m x bi, UNSPARSIFY collapses along the first dimension to give 
%     a full matrix which is the concatenation of Bi's along the second 
%     dimension, ie  X = [B1,B2,...,BN] and b is the vector of block sizes 
%     [b1,b2,...bN]. m must be the same for all blocks.
%       
%       UNSPARSIFY is the inverse function of SPARSEBLOCK.
% 
% [X,b] = unsparsify( SpX , 'transpose' ) transposes SpX and X so that SpX is
%     collapsed across columns rather than rows. This is the inverse
%     function of sparseblock( X, b, 'transpose' ).
% 
%     See also SPARSEBLOCK.

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

if isempty(SpX)
    error('Matrix must be non-empty')
end
if transpose
    SpX = SpX';
end

if nargout > 1
    [X,Ir,Jc] = unsparsifymex( SpX ); %unsparsifymex is a mex function

    npar = unique(diff(Jc));

    b = diff(find([1;diff(Ir(1:npar:end))>0;1]));
else
    X = unsparsifymex( SpX );
end

if transpose
    X = X';
end
