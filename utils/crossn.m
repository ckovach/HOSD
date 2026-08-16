
function tab = crossn(N,full)

% tab = crossn(N)
%
% Creates a matrix whose rows contain all combinations of integers
% 1:N(1),1:N(2),....
%
% For example, if N is a 2 X 1 vector of the lengths of 2 vectors, x1 and x2,
% then crossn creates a  N(1)*N(2) X 2 matrix which contains every pair of
% indices into x1 and x2.
%
%
% tab = crossn(N,0)
%
% Returns only combinations that are unique under permutation.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C. Kovach 2010


if nargin < 2 || isempty(full)
    full = true;
end


if length(N) == 1    
    tab = (1:N)';
    return
end
    
cn = crossn(N(2:end),full);
    
tab = cat(2,kron((1:N(1))',ones(size(cn,1),1)), repmat(cn,N(1),1));    
    
if ~full %Discard permutations of the same vector
    
%    [sqi,sqj]  = meshgrid(1:N(1),1:size(cn,1)); 
%    
%    keep = sqi <= sqj | sqi > size(cn,1);
   
   cols = 1:size(tab,2);
   
   [mn,mni] = min(max(tab));        

   keep = all(repmat(tab(:,mni),1,length(cols)-1) - tab(:,setdiff(cols,mni)) <=0,2);  
   
   tab = tab(keep,:);

end