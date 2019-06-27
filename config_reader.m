function config = config_reader(filename)
if nargin == 0
    filename = 'rnaseq_config.ini';
end

config = struct;
section = '';
temp = [];

fid = fopen(filename,'rt');

while 1
    line = fgetl(fid);
    if line == -1
        break;
    end
    if isempty(line)
        continue
    end
    
    if line(1) == '#'
        continue
    end
    if line(1) == '['
        if ~isempty(temp)
            % push the temp onto the main config, skip this the very first
            % time
            config = setfield(config,section,temp);
        end
        ind0 = strfind(line,'[');
        ind1 = strfind(line,']');
        section = line(ind0+1:ind1-1);   
        temp = struct;
        continue
    end
    
    
    % otherwise we just get the variables and assign the values
    parts = strsplit(line,'=');         
    temp = setfield(temp,strtrim(parts{1}),strtrim(parts{2}));
    
    
end

fclose(fid);
config = setfield(config,section,temp);