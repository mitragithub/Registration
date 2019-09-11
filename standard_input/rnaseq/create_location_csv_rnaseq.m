function create_location_csv_rnaseq(directory, dx, dy, thin_dz, thick_dz, thick_sections, outdir)
% because the sorting key may change, this should ONLY be used carefully
% also this file hard codes 15 micron pixels
if nargin < 7
    outdir = directory; % same input and output
end

files = dir([directory '*.tif']);
% we have to sort the files
key = cellfun(@(x)str2double(x(end-7:end-4)),{files.name},'uniformoutput',1);
% this sorting key may change
[sorted,permutation] = sort(key);
files = files(permutation);

count = 0;
% get the section id (last 4 numbers)
for f = 1 : length(files)
    % separate by underscores and get the last section
    [~,fname,~] = fileparts(files(f).name);
    ind = find(fname == '_',1,'last');
    count = count + 1;
    numbers(count) = str2num(fname(ind+1:end));
end

% get z locations.  Note that there may be missing files that correspond to
% some of these locations
n = max(numbers);
isthick = 0;
for i = 1 : n
    if i == 1
        % z will be the location of the start of the slice
        % c will be the location of the center of the slice
        z(1) = 0;         
        dz(1) = thin_dz;
        c(1) = thin_dz/2;
    else
        if any(i==thick_sections)
            z(i) = z(i-1) + thick_dz;    
            isthick(i) = 1;
            dz(i) = thick_dz;
        else
            z(i) = z(i-1) + thin_dz;
            isthick(i) = 0;
            dz(i) = thin_dz;
        end
        c(i) = c(i-1) + dz(i-1)/2 + dz(i)/2;
        
    end
end

% offset
c0 = (c(1) + c(end) )/2;
c = c - c0;
z = z - c0;


% now we just want my slices
zs = z(numbers);
cs = c(numbers);
isthicks = isthick(numbers);

% now we write it out
fid = fopen([outdir 'geometry.csv'],'wt');
% first print the fields
fprintf(fid,'filename, nx, ny, nz, dx, dy, dz, x0, y0, z0\n');
for i = 1 : length(files)
    fprintf(fid,'%s, ', files(i).name);
    info = imfinfo([directory files(i).name]);
    fprintf(fid,'%d, ', info.Width);
    fprintf(fid,'%d, ', info.Height);
    fprintf(fid,'%d, ', 1);
    fprintf(fid,'%f, ', dx);
    fprintf(fid,'%f, ', dy);
    if isthicks(i)
        fprintf(fid,'%f, ', thick_dz);
    else
        fprintf(fid,'%f, ', thin_dz);
    end
    
    % for offset we will subtract the mean
    x = (0 : info.Width-1)*dx;
    x = x - mean(x);
    fprintf(fid,'%f, ', x(1));
    y = (0 : info.Height-1)*dy;
    y = y - mean(y);
    fprintf(fid,'%f, ', y(1));
    % now for the z offset 
    fprintf(fid,'%f, ', cs(i));
    fprintf(fid,'\n');
end

% close the file
fclose(fid);