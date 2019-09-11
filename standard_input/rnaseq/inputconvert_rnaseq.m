%% 1. Select I/O directories
inputdir=input('Please specify the input directory (default to current): ', 's');
if isempty(inputdir)
    inputdir=pwd;
else
    if ~exist(inputdir,'dir')
    disp('Specified directory does not exsit! Default to the current directory.')
        inputdir=pwd;
    end
end
cd(inputdir)
outdir=input('Please specify the output directory (default to input directory): ', 's');
if isempty(outdir)
    outdir=inputdir;
else
    if ~exist(outdir,'dir')
     disp('Specified directory does not exsit! Default to the input directory.')
        outdir=inputdir;
    end
end
%% 2. Assemble geometry.csv
thickzinfofile=[inputdir,'/thickzinfo.csv'];
if ~exist('thickzinfofile','file')
    thickzinfofile=input('Please manually identify the csv file with thick sections information: ','s');
end
thickz=csvread(thickzinfofile); 
thick_dz=thickz(:,1);
thick_sections=thickz(:,2);
%
dx=.46*32;
dy=.46*32;
disp('Please verify the default values for x-y resolution of the input image: ')
disp(dx)
disp(dy)
resaccept=input('Confirm? (y/n) ','s');
if resaccept~='y'
    disp('Please specify the x-y resolution of the input image: ')
    dx=input('dx = (µm) ');
    dy=input('dy = (µm) ');
end
% 
thin_dz=20;
disp('Please verify the default values for z resolution of the input image: ')
disp(thin_dz)
resaccept=input('Confirm? (y/n) ','s');
if resaccept~='y'
    disp('Please specify the z resolution of the input image: ')
    thin_dz=input('dz = (µm) ');
end
create_location_csv_rnaseq(inputdir, dx, dy, thin_dz, thick_dz, thick_sections, outdir)