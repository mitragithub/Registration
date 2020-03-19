% write an image in vtk format
% use "simple legacy format"
% https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
% let my array I be of size
% ny x nx x nz x {1,3} x n_datasets
% where it has 1 component for scalars, or 3 components for vectors
%

function write_vtk_image(x,y,z,I,filename,title,names)
if nargin < 7
    names = {};
    for i = 1 : size(I,5)
        names{end+1} = ['dataset_' num2str(i-1)];
    end
end
if nargin < 6    
    title = 'no title specified';
    warning(['No title specified, using "' title '" as title'])
end

if length(names) ~= size(I,5)
    error(['Number of dataset names ' num2str(length(names)) ' does not match image size ' num2str(size(I,4)) '.'])
end

[directory, file, extension] = fileparts(filename);
if ~strcmp(extension, '.vtk')
    warning(['Filename does not have vtk extension, using provided extension ' extension ' anyway'])    
end

if size(I,4) == 1
    TYPE = 'SCALARS';
elseif size(I,4) == 3
    TYPE = 'VECTORS';
else
    error(['Image should be scalar or 3 component vector but has ' num2str(size(I,4)) ' components.'])
end


% open the vtk file
fid = fopen(filename,'wb');

% header line
fprintf(fid,'# vtk DataFile Version 2.0\n');

% title line
fprintf(fid,'%s\n',title);

% specify binary
fprintf(fid,'BINARY\n');

% type of dataset is structured points
fprintf(fid,'DATASET STRUCTURED_POINTS\n');

% dimension
% we use x,y,z as compared to matlabs row,col slice
dimension = [size(I,2),size(I,1),size(I,3)];
fprintf(fid,'DIMENSIONS %d %d %d\n',dimension);


% origin
origin = [min(x), min(y), min(z)];
fprintf(fid,'ORIGIN %f %f %f\n',origin);

% spacing
spacing = [x(2)-x(1), y(2)-y(1), z(2)-z(1)];
fprintf(fid,'SPACING %f %f %f\n',spacing);


% the data
num = size(I,1)*size(I,2)*size(I,3);
fprintf(fid,'POINT_DATA %d\n',num);



% loop through each image channel
% dataType is one of the types bit, unsigned_char, char, unsigned_short, short, unsigned_int, int,
% unsigned_long, long, float, or double.
type = class(I);
typestr = 'double';
if strcmp(type,'double')
    typestr = 'double';
elseif strcmp(type,'single')
    typestr = 'float';
elseif strcmp(type,'uint8');
    typestr = 'unsigned_char';
elseif strcmp(type,'uint16');
    typestr = 'unsigned_short';
elseif strcmp(type,'uint32');
    typestr = 'unsigned_int';
end
for i = 1 : size(I,5)
    if strcmp(TYPE,'SCALARS')
        fprintf(fid,'SCALARS %s %s %d\n', names{i}, typestr, 1); % note I'm using 1 component for each dataset
        fprintf(fid,'LOOKUP_TABLE default\n');
        % write this channel
        C = I(:,:,:,1,i);
        % make x fastest varying
        C = permute(C,[2,1,3,4]);
        fwrite(fid,C(:),type); % use matlab type
    elseif strcmp(TYPE, 'VECTORS')
        fprintf(fid,'VECTORS %s %s\n', names{i}, typestr);
        C = I(:,:,:,:,i);
        % make x fastest
        C = permute(C,[2,1,3,4]);
        % put vectors first
        C = permute(C,[4,1,2,3]);
        fwrite(fid,C(:),type); % use matlab type
    end
    
end


% close the vtk file
fclose(fid);