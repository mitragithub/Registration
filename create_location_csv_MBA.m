function create_location_csv_MBA(directory, dx, dy, dz, outdir,ext)
% because the sorting key may change, this should ONLY be used carefully
% also this file hard codes 15 micron pixels
% note it is actually 0.46*32 = 14.72
% width is 20, but distance from nissl to nissl is 40
% e.g. run like this
% create_location_csv_MBA('/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/PTM830/', 14.72, 14.72, 20)
% here we will identify the nissl images, and the "other" images

% keyboard
if nargin < 6
    ext = 'tif';
end
if nargin < 5
    outdir = directory; % same input and output
end
disp(['input directory ' directory])

% keyboard
files_ = dir([directory '*-N*.' ext]); % get the nissl only
files = {files_.name};
filesAll_ = dir([directory '*.' ext]); % get all files
filesAll = {filesAll_.name};
% remove nissl and keep what's left
filesF = setdiff(filesAll,files);

% the are numbered sequentially, take a half slice offset for nissl versus
% fluoro
% get the four digit number at the end for sorting
key = cellfun(@(x)str2double(x(end-7:end-4)),files,'uniformoutput',1);
keyF = cellfun(@(x)str2double(x(end-7:end-4)),filesF,'uniformoutput',1)+0.5;

files = [files,filesF];
key = [key,keyF];

[~,permutation] = sort(key);
files = files(permutation);
key = key(permutation);
z = key*dz*2;
z = z - (z(1) + z(end))/2;


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