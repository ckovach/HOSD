function out = cleangpu(me,in)

%Gather all fields that contain gpuArray objects
if nargin < 2
    in = me;
end

out = in;
if isstruct(in)
    fn = fieldnames(in);
else
    metac = metaclass(in);
    props = metac.PropertyList;
    props = props(~[props.Dependent]); %Exclude dependent properties that are actually get and set methods for some other property
    fn = setdiff({props.Name},{'Handle'});
end

if ~isempty(fn)
    for kk =1 :length(in)
        for k = 1:length(fn)
            try
                out.(fn{k})=me.cleangpu(in.(fn{k}));
            catch 
            end
        end
    end
else
    if iscell(in)
        for k = 1:length(in)
            out{k} = gather(in{k});
        end
    else
        out = gather(in);
    end
end

