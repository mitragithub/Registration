function [name,acronym] = allen_ontology_from_id(id);

% take id as an input, and search the file, and return name and acronym
if ~exist('allen_mouse_ontology.json','file')
    websave('allen_mouse_ontology.json','http://api.brain-map.org/api/v2/structure_graph_download/1.json');
end

fid = fopen('allen_mouse_ontology.json','rt');
while 1
    line = fgetl(fid);
    if line == -1
        break
    end
    
    if contains(line,'"id"')
        % get the id
        [T,R] = strtok(line,':');
        id_ = str2num(R(3:end-1));
        if id_ ~= id
            continue
        end
        % otherwise get new info
        line = fgetl(fid); % atlas id
        line = fgetl(fid); % ontology id
        line = fgetl(fid); % acronym
        [T,R] = strtok(line,':');
        acronym = R(4:end-2);
        line = fgetl(fid); % name
        [T,R] = strtok(line,':');
        name = R(4:end-2);
        break;
    end
end
fclose(fid);