function savevtksurface(filename,v,f,pointData,info)

if nargin < 4
    pointData = [];
end

% open the file
fid = fopen(filename,'wt');


% print the following lines
fprintf(fid,'# vtk DataFile Version 3.0\n');
fprintf(fid,'Surface Data\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET POLYDATA\n');
fprintf(fid,'\n');

% now we'll start writing points
nPoints = size(v,1);
fprintf(fid,'POINTS %d float\n',nPoints);
fprintf(fid,'%f %f %f\n',v');

% then we write the polygons
nPoly = size(f,1);
nEl = nPoly*4;
fprintf(fid,'POLYGONS %d %d\n',nPoly,nEl);
fprintf(fid,'3 %d %d %d\n',f'-1);

% now the point data
if ~isempty(pointData)
    fprintf(fid,'POINT_DATA %d\n',nPoints);
    nData = length(pointData);
    for i = 1 : nData
        if size(pointData(i).data,2) == 3
            fprintf(fid,'VECTORS %s float\n',pointData(i).name);
            fprintf(fid,' %f %f %f\n',pointData(i).data');% note the space
        else % scalars
            fprintf(fid,'SCALARS %s float 1\n',pointData(i).name);
            fprintf(fid,'LOOKUP_TABLE default\n');
            fprintf(fid,' %f\n',pointData(i).data);% note the space
        end
    end
end


% close the file
fclose(fid);