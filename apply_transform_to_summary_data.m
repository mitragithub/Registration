% the resolution of the image in registered space is 10x10x20
clear all;
close all;
fclose all;

% add functions for loading vtk files
% these addpath statements will depend on your local environment and will
% likely need to be changed
addpath /cis/home/dtward/Functions/vtk
addpath /cis/home/dtward/Functions/plotting


% filenames
% this should be changed for your local filesystem
geometry_file = 'geometry.csv';
image_file = 'PMD2159_cell_10.mat';


deformation_file = 'atlas_to_registered_displacement.vtk';
dz0 = 40; % expected slice sampling interval
dx0 = 10; 
output_name = 'data_registered_to_atlas.vtk';
output_title = 'data_registered_to_atlas';



% we will load the atlas to get a set of sample points for the output data
atlas_at_10_mu = '/cis/home/dtward/public_html/share/allen_vtk/average_template_10.vtk';

%%
% 1. load geometry file
fid = fopen(geometry_file,'rt');
count = 0;
z = [];
while 1
    line = fgetl(fid);
    if line == -1
        break
    end
    count = count + 1;
    
    if count == 1
        headings = line;
        continue
    end
    
    % get data, the only thing I need is the z0 coordinate
    % and only if it is fluoro
    data = strsplit(line,',');
    if isempty(strfind(data{1},'-F'))
        continue
    end
    z = [z,str2num(data{10})];
end
%%
% 2. load reconstructed volume, accounting for missing slices
vars = load(image_file); % this takes a really long time, probably there is some compression/decompression
J = vars.neurondensity;
clear vars;
%%
% check number of slices
% 
clear J_
nslice = size(J,3); % it is 255, same as the length of z
% find where ther are gaps
dz = diff(z);
inds = (z-z(1))/dz0+1;

figure;
imagesc(squeeze(sum(J,3)))

J_ = zeros(size(J,1),size(J,2),max(inds),size(J,4));
J_(:,:,inds,:) = J;
J = J_;
clear J_;

all_inds = 1 : size(J,3);
is_missing = zeros(size(all_inds));
for i = 1 : length(all_inds)
    is_missing(i) = ~sum(all_inds(i) == inds);
end
% fill in each missing slice with average of its neighbors,
% repeat several times in case there are multiple missing in a row
niter = 5;
for it = 1 : niter
    disp(['Filling in missing slices with neighbors ' num2str(it) ' of ' num2str(niter)])
    J_ = J;
    for i = 1 : size(J,3)
        if is_missing(i)
            if i == 1
                % set equal to next slice
                J(:,:,i,:) = J(:,:,i+1,:);
            elseif i == size(J,3)
                % set equal to previous slice
                J_(:,:,i,:) = J(:,:,i-1,:);
            else
                J_(:,:,i,:) = 0.5*(J(:,:,i-1,:) + J(:,:,i+1,:));
            end
        end
    end
    J = J_;
end
clear J_;

%%
% get the domain
% let z be defined by geometry
zJ = z(1) + (all_inds-1)*dz0;
% x is centered
xJ = (1 : size(J,2))*dx0; xJ = xJ - mean(xJ);
yJ = (1 : size(J,1))*dx0; yJ = yJ - mean(yJ);
%%
% 3. now we need to apply the deformation
% we want to go into atlas space at 10x10x20
[xI,yI,zI,deltaphi,title_,names] = read_vtk_image(deformation_file);
[XI,YI,ZI] = meshgrid(xI,yI,zI);
phi = deltaphi + cat(4,XI,YI,ZI);
% transform the image at 50 micron (not necessary)
% phiJ = zeros([size(phi,1),size(phi,2),size(phi,3),size(J,4)]);
% for c = 1 : size(J,4)
%     F = griddedInterpolant({yJ,xJ,zJ},J(:,:,:,c),'linear','none');
%     phiJ(:,:,:,c) = F(phi(:,:,:,2),phi(:,:,:,1),phi(:,:,:,3));
% end
% 

%%
% transform the image at 10x10x40
% I'll use the sampling grid of the atlas
% and downsampling the z
[xI_,yI_,zI_,I_,title_,names] = read_vtk_image(atlas_at_10_mu);
zI_ = zI_(1) : dz0 : zI_(end);
xI_ = xI_(1) : dx0 : xI_(end);
yI_ = yI_(1) : dx0 : yI_(end);
%%
[XI_,YI_,ZI_] = meshgrid(xI_,yI_,zI_);
% resample phi
phi_ = zeros([size(XI_),3]);
for c = 1 : 3
    F = griddedInterpolant({yI,xI,zI},phi(:,:,:,c),'linear','nearest');
    phi_(:,:,:,c) = F(YI_,XI_,ZI_);
end
% transform data
phiJ = zeros([size(phi_,1),size(phi_,2),size(phi_,3),size(J,4)]);
for c = 1 : size(J,4)
    F = griddedInterpolant({yJ,xJ,zJ},J(:,:,:,c),'linear','none');
    phiJ(:,:,:,c) = F(phi_(:,:,:,2),phi_(:,:,:,1),phi_(:,:,:,3));
end
phiJ(isnan(phiJ)) = 0;
%%
% write out in vtk format, here I will use single precision to save disk
write_vtk_image(xI_,yI_,zI_,single(phiJ),output_name,output_title,names)