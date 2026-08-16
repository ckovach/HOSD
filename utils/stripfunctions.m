
function R = stripfunctions(R)

% Removes all function handles from a structure array, descending through all fields. 
% One possible reason for wanting to do this is that saving a function handle to a .mat file 
% entails saving the entire workspace visible to the function, resulting in
% unnecessarily large files.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

%C. Kovach 2011


if isempty(R) || isnumeric(R)
    return
end


% if ~isstruct(R)
%     error('Input must be a structure or structure array, or else empty.')
% end
try
    fieldns = fieldnames(R);
catch
    return
end

for k = 1:length(fieldns)
    
    fn = fieldns{k};
    for i = 1:numel(R)
        try
            R(i).(fn);
        catch
            continue
        end
        if isa(R(i).(fn),'function_handle')  %If field contains a function handle, remove it           
            R(i).(fn) = [];  
%             R(i).(fn)(1) = [];          % Retain the field as a zero-length function handle
                                        % to avoid possible inconsistencies later.
        
        
        else %if isa( R(i).(fn) ,'struct')  %If field contains a structure run stripfunctions recursiveley
            try
                if ~isnumeric(R(i).(fn))
                    R(i).(fn) = stripfunctions( R(i).(fn) );          
                end
            catch 
                continue
            end
            
        end
        
    end 
    
end



