%% vtkimageoverlay.m
% This script generates image overlay between the raw image and the atlas
% segmentation
% call functions:
%   - read_vtk_image.m
%   - transparent_overlay.m (in mitragithub/DAP repository)
% Input:
%   - filename: string. file name only without path or extension
%   - imagedir: string. directory where the images are stored
%   - vtkdir: string. directory where the vtk outputs are stored
%   - vtkprefix: string. prefix of the vtk file
%   - dsrate: integer. downsample rate between the image and the vtk image,
%   typically vtk image is at equal or lower resolution.
% Output:
%   - imgpad: numeric matrix. origin image with padding added, if applicable.
%   - segimgup: resized segmentation map to match the image resolution
% Example:
% [jp2imgpad,segimgup]=vtkimageoverlay('M6328-F21--_2_0042',...
% '~/CSHLservers/mitragpu5/M25/marmosetRIKEN/NZ/m6328/m6328F/JP2/',...
% '~/Dropbox (Mitra Lab)/Data and Analysis/Marmoset/Connectivity/MarmosetBA/M6328/Native_Space_Outputs_v00/2D_nissl_native_space/',...
% 'atlas-seg_to_',64,1);
%
function [imgpad,segimgup]=vtkimageoverlay(filename,imagedir,vtkdir,vtkprefix,dsrate,ifvis)
% visualization flag
if nargin<6
    ifvis=0;
end
% read original image
imagefile=[imagedir,'/',filename,'.jp2'];
jp2img=imread(imagefile);
[xin,yin,~]=size(jp2img);
% read the vtk image
vtkfile=[vtkdir,'/',vtkprefix,'/',filename,'.vtk'];
[x,y,z,segimg,title,names,spacing,origin] = read_vtk_image(vtkfile); 
segimgup=imresize(segimg,dsrate,'nearest'); % upsample using nearest neighbor method
[xout,yout]=size(segimgup);
% pad the original image assuming it is placed in the center
imgpad=cast(zeros(xout,yout,3),'like',jp2img);
imgpad(xout/2-xin/2+1:xout/2+xin/2,yout/2-yin/2+1:yout/2+yin/2,:)=jp2img;
% visualization, if needed
if ifvis>0
    figure, imagesc(uint8(imgpad))
    transparent_overlay(gca,segimgup);
end