in_dir = 'DK39/DK39/';
out_dir = 'DK39/DK39RGB/';
if ~exist(out_dir,'dir')
    mkdir(out_dir)
end
files = dir([in_dir '*.tif']);

for i = 1 : length(files)
    I = imread([in_dir files(i).name]);
    
    
    I = repmat(I,[1,1,3]);
    imwrite(I,[out_dir files(i).name(1:end-4) '.png']);
    
    
end
