function [x,y,z,I,title_,names] = read_vtk_image(filename)
% note title_ has trailing underscore to not interfere with built in title
% function

% skip any white lines
if nargin < 1
    filename = 'test.vtk';
end
fid = fopen(filename,'rb');

% read the first line, it should just say 
% # vtk DataFile Version 2.0
while 1    
header = fgetl(fid);
if isempty(header) % skip any blank
    continue
end
disp(header)
break;
end


% now read the title
while 1    
title_ = fgetl(fid);
if isempty(title_) % skip any blank
    continue
end
disp(title_)
break;
end

% now read the type, it must be binary
while 1    
binary = fgetl(fid);
if isempty(binary) % skip any blank
    continue
end
disp(binary)
break;
end
if ~strcmp(binary,'BINARY')
    error(['Only support reading BINARY data, but datatype is ' binary])
end

% now dataset
while 1
dataset = fgetl(fid);
if isempty(dataset) % skip any blank
    continue
end
break;
end
[~,dataset] = strtok(dataset,' ');
dataset = dataset(2:end);
disp(dataset)
if ~strcmp(dataset,'STRUCTURED_POINTS')
    error(['Only support STRUCTURED_POINTS dataset, but dataset is ' dataset])
end


% now dimensions
while 1
dimensions = fgetl(fid);
if isempty(dimensions) % skip any blank
    continue
end
break;
end
[~,dimensions] = strtok(dimensions,' ');
dimensions = dimensions(2:end);
dimensions = sscanf(dimensions,'%d')';
disp('dimensions');disp(dimensions)
if length(dimensions) ~= 3
    error(['Only support 3D data, but dataset dimensions is ' num2str(dimensions)])
end

% now origin
while 1
origin = fgetl(fid);
if isempty(origin) % skip any blank
    continue
end
break;
end
[~,origin] = strtok(origin,' ');
origin = origin(2:end);
origin = sscanf(origin,'%f')';
disp('origin');disp(origin)
if length(origin) ~= 3
    error(['Only support 3D data, but dataset origin is ' num2str(origin)])
end

% now spacing
while 1
spacing = fgetl(fid);
if isempty(spacing) % skip any blank
    continue
end
break;
end
[~,spacing] = strtok(spacing,' ');
spacing = spacing(2:end);
spacing = sscanf(spacing,'%f')';
disp('spacing');disp(spacing)
if length(spacing) ~= 3
    error(['Only support 3D data, but dataset spacing is ' num2str(spacing)])
end

% create locations of data
x = (0:dimensions(1)-1)*spacing(1) + origin(1);
y = (0:dimensions(2)-1)*spacing(2) + origin(2);
z = (0:dimensions(3)-1)*spacing(3) + origin(3);


% now we will start point data
while 1
pointdata = fgetl(fid);
if isempty(pointdata) % skip any blank
    continue
end
break;
end
disp(pointdata)
POINT_DATA = 'POINT_DATA';
if ~strcmp(pointdata(1:length(POINT_DATA)), POINT_DATA)
    error(['Only support POINT_DATA, but found ' pointdata]);
end
[~,pointdata] = strtok(pointdata,' ');
pointdata = pointdata(2:end);
n_datapoints = sscanf(pointdata,'%d')';
disp(['Number of datapoints ' num2str(n_datapoints)])
if prod(dimensions) ~= n_datapoints
    error(['Product of dimensions should equal number of datapoints, but they are ' num2str(prod(dimensions)) ' and ' num2str(n_datapoints)])
end

% now we will begin looping through datasets
% I don't know how many datasets there are at the beginning
% so I'll use cells (better to add to end since don't need to be contiguous)
I = {};
names = {};
channels = {};
while 1
    line = fgetl(fid);
    % should say SCALARS or VECTORS
    if isempty(line)
        continue;
    end
    if line == -1
        break
    end
    [TYPE,R] = strtok(line, ' ');
    R = R(2:end);
    disp('image type')
    disp(TYPE)
    
    
    [names{end+1},R] = strtok(R,' ');
    R = R(2:end);
    disp('name')
    disp(names{end})
    
    [dtype,R] = strtok(R,' ');
    disp('dtype')
    disp(dtype)
    % we should do a conversoin
    dtype_matlab = 'uint8'; % default
    if strcmp(dtype,'unsigned_short')
        dtype_matlab = 'uint16';
    elseif strcmp(dtype,'unsigned_int')
        dtype_matlab = 'uint32';
    elseif strcmp(dtype,'float')
        dtype_matlab = 'single';
    elseif strcmp(dtype,'double')
        dtype_matlab = 'double';
    end
    
    if ~isempty(R)
        R = R(2:end);
        channels{end+1} = sscanf(R,'%d');
    else
        channels{end+1} = 1;
    end
    disp('channels')
    disp(channels{end})
    if channels{end} ~= 1
        error('only single channel images supported')
    end
    
    % now lookup table for scalars, I don't need it though
    if strcmp(TYPE,'SCALARS')
        line = fgetl(fid);
        [T,R] = strtok(line,' ');
        if ~strcmp(T,'LOOKUP_TABLE')
            error('after dataset type should be lookup table');
        end
    end
    
    % now read the data
    if strcmp(TYPE,'SCALARS')
        I{end+1} = fread(fid,n_datapoints,dtype_matlab);
        I{end} = reshape(I{end},dimensions);
        I{end} = permute(I{end},[2,1,3]); % switch x-y to row-col for matlab
    elseif strcmp(TYPE,'VECTORS')
%         keyboard
        I{end+1} = fread(fid,n_datapoints*3,dtype_matlab);
        I{end} = reshape(I{end},[3,dimensions]);
        I{end} = permute(I{end},[2,3,4,1]); % put components last
        I{end} = permute(I{end},[2,1,3,4]); % switch x-y to row-col
    end
end
fclose(fid);

nC = length(I);
I = reshape(I,[1,1,1,1,nC]);
I = cell2mat(I);

