function write_vtk_image(x,y,z,I,filename,title_,names)
% Write a vtk image
% Uses simple legacy vtk format which has human readable headers
% see here: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
% supports 1 (scalar) or 3 (vector) channel images, 
% supports multiple datasets in the same coordinate system 
% (e.g. multiple time points)
%
% inputs:
%     x,y,z                - spatial location of each voxel in the image,
%                            note x corresponds to the second image index 
%                            and y corresponds to the first image index
%     I                    - output image or vector field. Vector field
%                            uses the 4th index.  Multiple datasets (like 
%                            timeseries) uses the 5th index
%     filename             - a string containing the filename (full or 
%                            relative path) of vtk file
%     names                - cell array containing the name of each dataset
%                            Note we use the convention that if the name
%                            ends in (b) it is big endian.
%     title__              - title_ of the dataset (any user specified
%                            string).  Note it has an underscore to avoid a
%                            naming conflict with the built in function
%                            title_
%
% outputs:                 - writes a file with the specified filename
%      
%     
%
%
% Note: this file format does not specify "endian" in headers.  Most
% reader software (e.g. ITK snap, paraview) defaults to big endian, even 
% though many systems default to little endian.  This code will always
% write in big endian, and will put the symbol '(b)' in each dataset name.
% This is for backwards compatibility with previous code that did not
% specify endian.


if nargin < 7
    names = {};
    for i = 1 : size(I,5)
        names{end+1} = ['dataset_' num2str(i-1)];
    end
end
if nargin < 6    
    title_ = 'no title_ specified';
    warning(['No title_ specified, using "' title_ '" as title_'])
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

% title_ line
fprintf(fid,'%s\n',title_);

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
    if length(names{i}) >= 3 && strcmp(names{i}(end-2:end),'(b)')
        writename = names{i};
    else
        writename = [names{i} '(b)'];
    end
    if strcmp(TYPE,'SCALARS')        
        fprintf(fid,'SCALARS %s %s %d\n', writename, typestr, 1); % note I'm using 1 component for each dataset
        fprintf(fid,'LOOKUP_TABLE default\n');
        % write this channel
        C = I(:,:,:,1,i);
        % make x fastest varying
        C = permute(C,[2,1,3,4]);
        fwrite(fid,C(:),type,0,'b'); % use matlab type always write big
        % note default is little endian which does not work with itksnap or
        % paraviewcl
%         fwrite(fid,C(:),type,0,'b'); % use matlab type
    elseif strcmp(TYPE, 'VECTORS')
        fprintf(fid,'VECTORS %s %s\n', writename, typestr);
        C = I(:,:,:,:,i);
        % make x fastest
        C = permute(C,[2,1,3,4]);
        % put vectors first
        C = permute(C,[4,1,2,3]);
        fwrite(fid,C(:),type,0,'b'); % use matlab type
    end
    
end

% close the vtk file
fclose(fid);