function [v,f,pointData,info] = loadvtksurface(filename)

% filename = 'predictData/caudateContvsHighTCAbothSEG_FWER.vtk';

pointDataCount = 0;
pointData = struct();
info = [];

% open the file for reading
fid = fopen(filename);

% let's read through line by line
while 1
    line = fgetl(fid);
    if line == -1 
        break;
    end
    
    % get the first word on the line
    [T,R] = strtok(line);
    
    % if it is points
    if strcmpi(T,'POINTS')
        % how many points
        [T,R] = strtok(R);
        nPoints = str2num(T);
        
        data = fscanf(fid,'%f %f %f\n',[3 nPoints]);
        v = data';
        
    end
    
    % if it is polygons
    if strcmpi(T,'POLYGONS')
        % how many points
        [T,R] = strtok(R);
        nPolys = str2num(T);
        % how much data
        [T,R] = strtok(R);
        nData = str2num(T);
        % I only support triangulated surfaces
        if nData/nPolys ~= 4
            error('Not a triangulated surface')
        end
        data = fscanf(fid,'%d %d %d %d\n',[4 nPolys]);
        f = data(2:4,:)' + 1;
    end
    
    
    
    % pointData
    if strcmpi(T,'POINT_DATA')
        % once we start reading point data, we keep reading it

        % how much data (should be same as number of vertices)
        nPointData = str2num(R);
        if nPointData ~= nPoints
            error('point data is not same size as points')
        end
        
        
        
        
        while 1
            % get the next line (watch for double spaces)
            line = '';
            while isempty(line)
                line = fgetl(fid);
                if line == -1
                    break;
                end
            end
            if line == -1 
                break;
            end
            
            % get the first word
            [T,R] = strtok(line);
            if strcmpi(T,'SCALARS')
                
                
                pointDataCount = pointDataCount + 1;
                
                % get the second word
                [T,R] = strtok(R);
                pointData(pointDataCount).name = T;
                
                % don't care about the next word
                
                % get the next line (watch out for double spaces)
                line = '';
                while isempty(line)
                    line = fgetl(fid);
                end
                % I don't care about this line
                
                % now get the data
                data = [];
                while isempty(data)
                    data = fscanf(fid,'%f\n',nPointData);
                end
                pointData(pointDataCount).data = data;
            elseif strcmpi(T,'VECTORS')
                
 
                pointDataCount = pointDataCount + 1;
                
                % get the second word
                [T,R] = strtok(R);
                pointData(pointDataCount).name = T;
                
                % don't care about the next word
                
                % get the next line (watch out for double spaces)
                line = '';
                while isempty(line)
                    line = fgetl(fid);
                end
                % I don't care about this line
                
                % now get the data
                data = [];
                while isempty(data)
                    data = fscanf(fid,'%f %d %f\n',nPointData*3);
                end
                pointData(pointDataCount).data = reshape(data,3,[])';
            end

        end
    
    end
    
    
end














% close the file
fclose(fid);




info = '';