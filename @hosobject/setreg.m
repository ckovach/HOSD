 function me=setreg(me,xin)
 
% Set the regressor property
   if isa(xin,'regressor')
        if isempty(xin)
            return
        end
        R = xin;
        xin = R.value;
    else
        R = [];
    end
    if size(xin,1)> me(1).buffersize
        for k = 1:size(xin,2)
            X(:,k) = sum(me(1).chop_input(xin(:,k)),1);
        end
    elseif size(xin,1)==me(1).buffersize
        for k = 1:size(xin,3)
            X(:,k) = sum(xin(:,:,k));
        end
    else
        X = xin;
    end
    if isempty(R)
        if isempty(which('regressor'))
            return
        else
            R = regressor(X);
        end
    else
        R2 = regressor(X);
        R.value = R2.value;
        if isfield(R,'m')
            R.m=[];
        end
%                 R.m = R2.m;
    end
    for k = 1:length(me)
        me(k).regval = R;
    end
 end