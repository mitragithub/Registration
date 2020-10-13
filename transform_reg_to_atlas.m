% Feb 2020, converted into function from apply_transform_to_summary_data.m by Daniel Tward
% each 10um transformation takes 15-20min
% profile is available at https://www.dropbox.com/sh/8k4xn8ooaqjyrma/AACgzyW9bz4gakluhy-TBw-la?dl=0
function phiJ=transform_reg_to_atlas(INPUTvol,INPUTspacing,geometry_file,deformation_file,atlas_vtk,output_name)
% % the resolution of the image in registered space is 10x10x20
% clear all;
% close all;
% fclose all;
%
% % add functions for loading vtk files
% % these addpath statements will depend on your local environment and will
% % likely need to be changed
% addpath /cis/home/dtward/Functions/vtk
% addpath /cis/home/dtward/Functions/plotting


% filenames
% this should be changed for your local filesystem
% geometry_file = 'geometry.csv';
% image_file = 'PMD2159_cell_10.mat';
% deformation_file = 'atlas_to_registered_displacement.vtk';
dz0 = INPUTspacing(3); % expected slice sampling interval
dx0 = INPUTspacing(1);
% output_name = 'data_registered_to_atlas.vtk';
% output_title = 'data_registered_to_atlas';
tic;


% we will load the atlas to get a set of sample points for the output data
% atlas_vtk = '~/Dropbox (Mitra Lab)/Data and Analysis/Mouse/MouseBrainAtlases/AllenMouseBrainAtlas/Standardized/Average_Template/average_template_10.vtk';

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
fclose(fid);
%%
% 2. load reconstructed volume, accounting for missing slices
% vars = load(image_file); % this takes a really long time, probably there is some compression/decompression
% INPUTvol = vars.neurondensity;
% clear vars;
%%
% check number of slices
%
clear J_
nslice = size(INPUTvol,3); % it is 255, same as the length of z
% find where ther are gaps
dz = diff(z);
inds = (z-z(1))/dz0+1;
%
% figure;
% imagesc(squeeze(sum(INPUTvol,3)))

J_ = zeros(size(INPUTvol,1),size(INPUTvol,2),max(inds),size(INPUTvol,4));
J_(:,:,inds,:) = INPUTvol;
INPUTvol = J_;
clear J_;

all_inds = 1 : size(INPUTvol,3);
is_missing = zeros(size(all_inds));
for i = 1 : length(all_inds)
    is_missing(i) = ~sum(all_inds(i) == inds);
end
% fill in each missing slice with average of its neighbors,
% repeat several times in case there are multiple missing in a row
niter = 5;
for it = 1 : niter
    disp(['Iteratively filling in missing slices with neighbors (' num2str(it) ' of ' num2str(niter),')'])
    J_ = INPUTvol;
    for i = 1 : size(INPUTvol,3)
        if is_missing(i)
            if i == 1
                % set equal to next slice
                INPUTvol(:,:,i,:) = INPUTvol(:,:,i+1,:);
            elseif i == size(INPUTvol,3)
                % set equal to previous slice
                J_(:,:,i,:) = INPUTvol(:,:,i-1,:);
            else
                J_(:,:,i,:) = 0.5*(INPUTvol(:,:,i-1,:) + INPUTvol(:,:,i+1,:));
            end
        end
    end
    INPUTvol = J_;
end
clear J_;

%%
% get the domain
% let z be defined by geometry
zJ = z(1) + (all_inds-1)*dz0;
% x is centered
xJ = (1 : size(INPUTvol,2))*dx0; xJ = xJ - mean(xJ);
yJ = (1 : size(INPUTvol,1))*dx0; yJ = yJ - mean(yJ);
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
[xI_,yI_,zI_,I_,title_,names] = read_vtk_image(atlas_vtk);
% zI_ = zI_(1) : dz0 : zI_(end);
% dxtarget=25;
% zI_ = zI_(1) : dxtarget : zI_(end); % convert into isotropic
% xI_ = xI_(1) : dxtarget : xI_(end);
% yI_ = yI_(1) : dxtarget : yI_(end);
%% this is the time consuming part
[XI_,YI_,ZI_] = meshgrid(xI_,yI_,zI_); % target space
% resample phi
phi_=cell(3,1); % each color channel is saved in one cell
for c = 1 : 3
    phi_{c} = zeros(size(XI_));
    F = griddedInterpolant({yI,xI,zI},phi(:,:,:,c),'linear','nearest');
    phi_{c} = F(YI_,XI_,ZI_); % evaluate phi at each point in target space
end
% transform data
phiJ=cell(3,1);
for c = 1 : size(INPUTvol,4)
    tic;
    phiJ{c} = zeros(size(phi_{1}));
    if sum(sum(sum(INPUTvol(:,:,:,c))))>0
    F = griddedInterpolant({yJ,xJ,zJ},INPUTvol(:,:,:,c),'linear','nearest');
    phiJ{c} = F(phi_{2},phi_{1},phi_{3});
    phiJ{c}(isnan(phiJ{c})) = 0;
    end
    toc;
end
phiJ1=cat(4,phiJ{1},phiJ{2},phiJ{3});
phiJ=phiJ1;
clear phiJ1
%%
if nargin>5 % save
    disp(['Saving ',output_name,'...'])
    neurondensityatlas=phiJ;
    if exist([output_name,'.mat'],'file')
        save([output_name,'.mat'],'neurondensityatlas','-append')
    else
        save([output_name,'.mat'],'neurondensityatlas','-v7.3')
    end
    % write out in vtk format, here I will use single precision to save disk
    %     write_vtk_image(xI_,yI_,zI_,single(phiJ),[output_name,'.vtk'],output_title,names)
    disp('Saved.')
    toc;
end