function out = connectmap(mask)

if size(mask,1) == 1
    mask = true(mask);
end

sz = size(mask);

cind = arrayfun(@(x)1:x,sz,'uniformoutput',false);

IJK = fliplr(crossn(fliplr(sz)));
if ~isscalar(mask)
    IJK = IJK(mask(:),:);
end

n = prod(sz);
m = size(IJK,1);


nbr = [-1 0 1];
cr = crossn(length(nbr)*ones(size(sz)));

cr(all(nbr(cr)==0,2),:) = [];

nindx = (IJK(:,1)+nbr(cr(:,1)));
nindx(nindx<1 | nindx>sz(1)) = nan;

for k = 2:length(sz)
    
    didx = (IJK(:,k)+nbr(cr(:,k)));
    didx(didx<1 | didx >sz(k)) = nan;
    nindx = (didx-1).*prod(sz(1:k-1)) + nindx;
end

bigas = (1:m)'+ m*(nindx-1);
Adj = spalloc(m,n,n*3^length(sz));
Adj(bigas(~isnan(bigas))) = 1;
Adj = Adj(:,mask(:));
out.Adj = Adj;

if ~all(mask(:))
    gr = graph(Adj);
    
    clusters = zeros(size(mask));
    clusters(mask) = conncomp(gr);
    out.clusters = clusters;
    out.graph = gr;
    % %Adj = mask(:).*Adj.*mask(:)';
    % 
    % varargout{1} = Adj;
    % 
    % %if nargout > 1
    %     nsing = sum(Adj)>0;
    % 
    %  fmask = find(mask(:));
    %     fmask = fmask(nsing);
    %     Adj = Adj(nsing,nsing);
    %     Dg = diag(sum(Adj));
    %     Dgi = Dg; 
    %     Dgi(find(Dg)) = 1./sqrt(Dg(find(Dg)));
    %     %Lap = Dg-Adj;
    %     %Lnorm = speye(size(Adj))-Dgi*Adj*Dgi;
    %     [u,l] = svd(full(Lnorm));
    %     nclust = size(u,1)-rank(l);
    %     getcl = size(u,1)-(nclust:-1:1)-nclust+1;
    %     varargout{2} = Lorm;
    % %end
end