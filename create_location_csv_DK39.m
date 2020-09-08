function create_location_csv_DK39(directory, dx, dy, dz, outdir, ext)
% use metadata and xlsx file to build


% keyboard
if nargin < 6
    ext = 'tif';
end
if nargin < 5
    outdir = directory; % same input and output
end
disp(['input directory ' directory])


fid = fopen('DK39/DK39.sections.csv','rt');
count = 0;
while 1
    line = fgetl(fid);
    if line == -1
        break
    end
    out = strsplit(line,',');
    if strcmp(out{1},'File')
        continue
    end
    
    count = count + 1;
    
    name{count} = out{1};
    num(count) = str2num(out{2});
end
fclose(fid);



% get a z location for each slice
z = num*dz;
z = z - (z(1) + z(end))/2;

files_ = dir([directory '*' ext]); 
% put them in the same order as the name
for i = 1 : length(name)
    ind = cellfun(@(x) ~isempty(x),strfind({files_.name},name{i}(1:end-4)));

    files{i} = files_(ind).name;
end



% now we write it out
fid = fopen([outdir 'geometry.csv'],'wt');
% first print the fields
fprintf(fid,'filename, nx, ny, nz, dx, dy, dz, x0, y0, z0\n');
for i = 1 : length(files)
    fprintf(fid,'%s, ', files{i});
    info = imfinfo([directory files{i}]);
    fprintf(fid,'%d, ', info.Width);
    fprintf(fid,'%d, ', info.Height);
    fprintf(fid,'%d, ', 1);
    fprintf(fid,'%f, ', dx);
    fprintf(fid,'%f, ', dy);
    fprintf(fid,'%f, ', dz);
    
    % for offset we will subtract the mean
    x = (0 : info.Width-1)*dx;
    x = x - mean(x);
    fprintf(fid,'%f, ', x(1));
    y = (0 : info.Height-1)*dy;
    y = y - mean(y);
    fprintf(fid,'%f, ', y(1));
    % now for the z offset 
    fprintf(fid,'%f, ', z(i));
    fprintf(fid,'\n');
end

% close the file
fclose(fid);