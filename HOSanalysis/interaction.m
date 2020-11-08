function RX = interaction(varargin)

%
%   Creates a regressor of interaction terms between R1 and R2. 
% 
% Usage:  RX = interaction(R1,R2)
% 
%         RX = interaction(R,'intxnord',intxnord)
%
%   Creates interactions among all the regressors within the regressor
%   array, R, up to order intxnord.By default intxnord = 2.
%
% function RX = interaction(R)
%
%   Creates 2nd order interactions among all the regressors within the regressor
%   array, R.
%
% See also PPREG

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/GazeReader/0.1/interaction.m $
% $Revision: 151 $
% $Date: 2014-03-31 11:10:53 -0500 (Mon, 31 Mar 2014) $
% $Author: ckovach $
% ------------------------------------------------

% C. Kovach 2007 - 20019
intxnord = 2;

if isstruct(varargin{1}) || isa(varargin{1},'regressor')
    R1 = varargin{1};
    varargin(1)=[];
end


if ~isempty(varargin) && (isstruct(varargin{1})|| isa(varargin{1},'regressor'))
    R2 = varargin{1};
    varargin(1)=[];
elseif ~isempty(varargin) && isnumeric(varargin{1})
%     R2 = varargin{1};
    intxnord = varargin{1};
else
    R2 = 0;
end



i=1;
codeincr = 0;
label = [];
maxord = Inf;
while i <= length(varargin)
    
   switch lower(varargin{i})
        case 'intxnord' %order of interaction 
            intxnord = varargin{i+1};
            i = i+1;            
        case 'codeincr'
            codeincr = varargin{i+1};
            i = i+1;            
        case 'maxord'  %Returns only interactions of order less than or equal to this term, as determined from the info.intxnord field.
            maxord = varargin{i+1};
            i = i+1;            
        case 'label'
            label = varargin{i+1};
            i = i+1;            
        otherwise
            error('%s is not a valid keyword.',varargin{i})
    end
    i = i+1;
    
end

if intxnord == 1 %interaction order 1 implies no terms other than the ones already in the input
    RX = R1([]);
    return
end


if nargin > 1 && (isstruct(R2)||isa(R2,'regressor'))
    %For two regressor arrays 2nd order interactions between groups
    indx = 1;
    RX = R2([]);% RX(1) = [];
    for i = 1:length(R1)
        for j = 1:length(R2)

            if any(R1(i).noptions(:) ~= R2(j).noptions(:))
                error('Blocks do not appear to be identical for all regressors')
            end
            int = interaction2(R1(i), R2(j),maxord);
            if ~isempty(int)
                RX(indx) = int;
                RX(indx).code = codeincr + indx;
                RX(indx).codevec(1:RX(indx).Npar) = codeincr + indx;
                if ~isempty(label)
                    RX(indx).label = sprintf('%s %i',label,indx);
                end
                indx = indx+1;
            end
        end
    end
    
    return
    
elseif length(R1) == 2
    
    RX = interaction(R1(1),R1(2),'codeincr',codeincr);
    RX.code = codeincr + 1;
    if ~isempty(label)
         RX.label = label;
    end
    return
    
elseif length(R1) == 1 && (nargin == 1  || isempty(R2) )
    
    RX = R1;
    return
    
end

%For one array, all interactions among members are computed up to specified
%interaction order.



indx = 1;
RX=R1([]);
for intxn = 2:intxnord
    B = chooseperm(length(R1),intxn);
    
    if isempty(B)
        break
    end
    for i = 1:size(B,1)
        RXsub =interaction(R1(B(i,1)), interaction( R1(B(i,2:end)) ) );
        for k = 1:length(RXsub)
            RXsub(k).code = codeincr + indx;
            RXsub(k).codevec(1:RXsub(k).Npar) = codeincr + indx;
           indx = indx+1;
        end
        RX = [RX,RXsub];
    end
    
end

% if nargin < 3
%     ixnlev = 2;
% end

try
    fid = fopen([mfilename,'.m'],'r');
    R.info.COMMAND =  fread(fid);
    fclose(fid);
catch %#ok<CTCH>
    warning('Failed to record the script used to generate regressor %s poly',R.label);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RX = interaction2(R1,R2,maxord)


rxind = 1;


RX = R1;
info = RX.info;

% RX = RX([]);
    
for i = 1:length(R1)
    if strcmp(R1.info.form,'sparse')
        Vi = unsparsify(R1.value,'transpose');
    else
        Vi = R1.value;
    end        
    
    for j = 1:length(R2)
        
        if R1(i).info.intxnord + R2(j).info.intxnord > maxord  %only return regressors of order less than maxord
            continue
        end
        RX(rxind).info = info;
        if strcmp(R2.info.form,'sparse')
            Vj = unsparsify(R2.value,'transpose');
        else
            Vj = R2.value;
        end        
    
        
%         RX(rxind).value = sparseblock(repmat( Vi,1,size( Vj,2)).*kron( Vj ,ones(1,size( Vi,2))), R1(i).noptions,'transpose');     
%         RX(rxind).levmat = cat(1,repmat( R1.levmat,1,size( Vj,2)), kron( R2.levmat ,ones(1,size( Vi,2))));
%         RX(rxind).factmat = cat(1,repmat( R1.factmat,1,R2.Npar), kron( R2.factmat + 1,ones(1,R1.Npar)));

        if strcmp(R2.info.form,'sparse') ||strcmp(R1.info.form,'sparse')  
            RX(rxind).value = sparseblock(kron( Vi,ones(1,size( Vj,2))).*repmat( Vj ,1,size( Vi,2)), R1(i).noptions,'transpose');     
        elseif ~isempty(Vi) && ~isempty(Vj)
            RX(rxind).value = kron( Vi,ones(1,size( Vj,2))).*repmat( Vj ,1,size( Vi,2));     
        else
            RX(rxind).value =[];
            RX(rxind).fun = makemefun([R1,R2]);
%             continue
        end
        RX(rxind).levmat = cat(1,kron( R1.levmat,ones(1,R2.Npar)), repmat( R2.levmat ,1,R1.Npar));
        RX(rxind).factmat = cat(1,kron( R1.factmat,ones(1,R2.Npar)), repmat( R2.factmat,1,R1.Npar));
        RX(rxind).normconst = R1(i).normconst*R2(j).normconst; 
        RX(rxind).label = sprintf('%s * %s',R1(i).label,R2(j).label) ;
        RX(rxind).info.label = sprintf('%s * %s',R1(i).label,R2(j).label) ;
        RX(rxind).info.intxnord = R1(i).info.intxnord + R2(j).info.intxnord;
        RX(rxind).info.parent = cat(1,R1(i).info,R2(j).info);
%         RX(rxind).info.hashcode = R1(i).info.hashcode./2 + R2(j).info.hashcode;
        if strcmp(R2.info.form,'sparse') ||strcmp(R1.info.form,'sparse')  
            RX(rxind).info.form = 'sparse';
        else
            RX(rxind).info.form = R1.info.form ;
        end
        RX(rxind).noptions= R1(i).noptions;
        RX(rxind).Npar = R1(i).Npar*R2(j).Npar;
        RX(rxind).fixed = zeros(1,RX(rxind).Npar);
%         [RX(rxind).function,RX(rxind).deriv,terms] = makemefun([R1,R2]);
%         RX.info.functionInputCodes = terms;
        
        rxind = rxind+1;
    end
end

%%%%%%%%%%%%%%%%%

function [mkfun,mkdfun,terms] = makemefun(Rs)

rcodes = [Rs.code]; %#ok<NASGU>
rfunctions = arrayfun(@(R)@(varargin)R.makefun(varargin{:}),Rs,'uniformoutput',false); %#ok<NASGU>
factmats = {Rs.factmat};
% dfunctions = {Rs.deriv}; %#ok<NASGU>

npar1 = Rs(1).Npar; %#ok<NASGU>
npar2 = Rs(2).Npar; %#ok<NASGU>


terms= zeros(2,0);
for i = 1:length(Rs)
    if ~isempty(factmats{i})
        [unq,~,polyterm] = unique(factmats{i}(:,1),'stable');
    else
        unq = [];
    end
     rinputs{i} = []; %#ok<*AGROW>
     for j = 1:length(unq)
         for k = 1:sum(polyterm==j)
             
             input = find(terms(1,:) == unq(j) & terms(2,:) == k);
             if isempty(input)
                 terms(:,end+1) = [unq(j),k];
                 rinputs{i}(end+1) = size(terms,2);
             else
                 rinputs{i}(end+1) = input;
             end
                 
         end
     end
end             

nargs = max([rinputs{:}]);

argnum = cellfun(@(n) sprintf('%i',n),num2cell(1:nargs),'uniformoutput',false);
args = strcat('X',argnum);

trim = @(a) a(1:end-1);
genarg = @(X)trim(sprintf('%s,',X{:}));
funs = cellfun(@(inp,i) sprintf('rfunctions{%i}(%s)',i,genarg(args(inp))),rinputs,num2cell(1:length(rinputs)),'uniformoutput',false); %#ok<NASGU>

dfuns = cellfun(@(inp,i) sprintf('dfunctions{%i}(%s,n([%s]))',i,genarg(args(inp)),genarg(argnum(inp))),rinputs,num2cell(1:length(rinputs)),'uniformoutput',false); %#ok<NASGU>

funstr = sprintf('@(%s) kron( rfunctions{1}(%s),ones(1,npar2) ).*repmat( rfunctions{2}(%s) , 1 , npar1 );',...
            genarg(args),...
            genarg(args(rinputs{1})),...
            genarg(args(rinputs{2})));

mkfun = eval(funstr);



try
    dfunstr = sprintf('@(%s,n) kron( dfunctions{1}(%s,n([%s])),ones(1,npar2) ).*repmat( dfunctions{2}(%s,n([%s])) , 1 , npar1 );',...
                genarg(args),...
                genarg(args(rinputs{1})),genarg(argnum(rinputs{1})),...
                genarg(args(rinputs{2})),genarg(argnum(rinputs{2})));

    mkdfun = eval(dfunstr);
catch %#ok<CTCH>

    mkdfun = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = chooseperm(N, k, item, rows)

% 
%    Returns indices for every unordered combination of k items from a
%    population of N using a recursive algorithm.


if nargin < 3
    rows = 1;
    item = 1;
end

B = [];

if item >= k
    
    B = (rows:N)';
    
else
    
    for strow = rows:N

       b = chooseperm( N, k, item+1, strow+1);
    
       if ~isempty(b)
           
           B = cat(1,B,cat(2,ones(size(b,1),1)*strow,b));
           
       else
           
           return
           
       end
    end    
end

