function transform_seg(seg_file, input_dir, output_dir, save_dir, res)

% seg file is the filename of the atlas you want to deform (vtk)
% input_dir is a directory of tif images you want to match, with a
% geometry csv file in the same directory
% output_dir is the directory that contains the outputs from registration
% res is the desired resolution for the data (e.g. microns per pixel)
%
% this function has no outputs, it just makes a picture.  You can save
% outputs by editing the bottom

% 

addpath Functions/vtk

addpath Functions/plotting
% if nargin == 0
%     % example inputs
%     seg_file = 'annotation_50.vtk';
%     input_dir = 'MD720/';
%     output_dir = 'MD720_OUT/';
%     res = 5; % e.g. 5 um per pixel
% end


%%
% load the atlas
[xI,yI,zI,I,title_,names] = read_vtk_image(seg_file);
rng(1);
colors = rand(256,3);
colors(1,:) = 1; % set background white

%% 
% get geometry of files from input dir
geometry_file = dir([input_dir '*.csv']);
fid = fopen([input_dir geometry_file(1).name],'rt');
line = fgetl(fid); % ignore the first line
% it should say
% filename, nx, ny, nz, dx, dy, dz, x0, y0, z0

csv_data = {};
count = 0;
while 1
    line = fgetl(fid);
    if line == -1
        break
    end
    count = count + 1;
    % process this line, splitting at commas
    csv_data(count,:) = strsplit(line,',');
    %     
end
fclose(fid);
files = csv_data(:,1);
nxJ = cellfun(@(x)str2num(x), csv_data(:,2:3));
dxJ = cellfun(@(x)str2num(x), csv_data(:,5:6));

x0J = cellfun(@(x)str2num(x), csv_data(:,8:9));
z0J = cellfun(@(x)str2num(x), csv_data(:,10));

zJ = z0J;
dzJ = cellfun(@(x) str2num(x), csv_data(:,7));

for f = 1 : length(files)
    xJ{f} = x0J(f,1) + (0:nxJ(f,1)-1)*dxJ(f,1);
    yJ{f} = x0J(f,2) + (0:nxJ(f,2)-1)*dxJ(f,2);
end


%%
% loop through the files
% for each file apply one transform to the atlas, and one to the target, so
% they meet in registered space

for f = 1 : length(files)
    [directory,fname,ext] = fileparts(files{f});

    % load the tif image and display
    J = imread([input_dir files{f}]);
    if strcmp(class(J),'uint8')
        J = double(J)/255.0;
    end
    J = double(J);
    
    % danfigure(1);
    % imagesc(xJ{f},yJ{f},J);
    % axis image;
    % title(['Input data ' strrep(fname,'_','\_')])
    
    %%
    % transform this image to registered space (2d transforms, identity in z)
    % this involves two steps
    % 1. upsample the transform to resolution res
    % 2. resample the data
    % 
    % load the transform
    try
        [xTJ,yTJ,zTJ,DeltaTJ,title_,names] = read_vtk_image([output_dir 'registered_to_input_displacement_' fname '.vtk']);
    catch
        disp(['Could not read ' fname])
        continue
    end
    
    % upsample the transform
    xTJup = xTJ(1) : res : xTJ(end);
    yTJup = yTJ(1) : res : yTJ(end);
    [XTJup,YTJup] = meshgrid(xTJup,yTJup);
    DeltaTJup = zeros(size(XTJup,1),size(XTJup,2),1,2);
    F = griddedInterpolant({yTJ,xTJ},DeltaTJ(:,:,1,1),'linear','nearest');
    DeltaTJup(:,:,1,1) = F(YTJup,XTJup);
    F = griddedInterpolant({yTJ,xTJ},DeltaTJ(:,:,1,2),'linear','nearest');
    DeltaTJup(:,:,1,2) = F(YTJup,XTJup);
    
    phiTJup = zeros(size(DeltaTJup));
    phiTJup(:,:,:,1) = DeltaTJup(:,:,:,1) + XTJup;
    phiTJup(:,:,:,2) = DeltaTJup(:,:,:,2) + YTJup;
    
    
    
    % apply transform to data    
    JT = zeros(size(phiTJup,1),size(phiTJup,2),size(J,3));
    for c = 1 : size(J,3)
        F = griddedInterpolant({yJ{f},xJ{f}},J(:,:,c),'linear','nearest');
        JT(:,:,c) = F(phiTJup(:,:,:,2),phiTJup(:,:,:,1));
    end
    
    % danfigure(2);
    % imagesc(xTJ,yTJ,JT)
    % axis image
    % title(['Registered data ' strrep(fname,'_','\_')])
    
    
    %%
    % transform the atlas
    % load the transform
    [xTI,yTI,zTI,DeltaTI,title_,names] = read_vtk_image([output_dir 'registered_to_atlas_displacement_' fname '.vtk']);
    
    % upsample it
    DeltaTIup = zeros(size(XTJup,1),size(XTJup,2),1,3);

    F = griddedInterpolant({yTJ,xTJ},DeltaTI(:,:,1,1),'linear','nearest');
    DeltaTIup(:,:,1,1) = F(YTJup,XTJup);
    F = griddedInterpolant({yTJ,xTJ},DeltaTI(:,:,1,2),'linear','nearest');
    DeltaTIup(:,:,1,2) = F(YTJup,XTJup);
    F = griddedInterpolant({yTJ,xTJ},DeltaTI(:,:,1,3),'linear','nearest');
    DeltaTIup(:,:,1,3) = F(YTJup,XTJup);
    
    phiTIup = zeros(size(DeltaTJup));
    phiTIup(:,:,:,1) = DeltaTIup(:,:,:,1) + XTJup;
    phiTIup(:,:,:,2) = DeltaTIup(:,:,:,2) + YTJup;
    phiTIup(:,:,:,3) = DeltaTIup(:,:,:,3) + zJ(f);

    % apply transform to data
    [XTI,YTI,ZTI] = meshgrid(xTI,yTI,zTI);
    
    
    F = griddedInterpolant({yI,xI,zI},I,'nearest','nearest'); % for annotation use nearest
    IT = F(phiTIup(:,:,:,2),phiTIup(:,:,:,1),phiTIup(:,:,:,3));
    save(strcat(save_dir, fname,'.mat'), 'IT')

    % for display
    % ITm = mod(IT,256);
    % IT_RGB = reshape([colors(ITm+1,1),colors(ITm+1,2),colors(ITm+1,3)], size(ITm,1),size(ITm,2),3);
    % danfigure(3);
    % imagesc(xTJ,yTJ,IT_RGB)
    % axis image
    % title(['Atlas transformed to ' strrep(fname,'_','\_')])
    
    %%
    % overlay them (assume that J is RGB in the range [0,1]
    % alpha = 0.25; % opacity of labels
    % Ishow = bsxfun(@plus, IT_RGB*alpha, JT*(1-alpha));
    % danfigure(4);
    % imagesc(xTJ,yTJ,Ishow)
    % axis image
    % title([fname])
    
    
%     %%
%     % show a deformation grid
%     levels = linspace(-10000,10000,500);
%     hold on;
%     contour(xTJup,yTJup,phiTIup(:,:,:,1),'color','r');
%     contour(xTJup,yTJup,phiTIup(:,:,:,2),'color','g');
% %     contour(xTJup,yTJup,phiTIup(:,:,:,3),'color','b'); % you may or may not want to look at slices
%     hold off;
    
%     %%
%     % you may want to save some of these figures
%     drawnow;

    
end
