function imgstack=vtk_image_batch(vtkdir,stackfilename,filenamepart,modality,idloc)
%% example: nisslstack=vtk_image_batch('2D_nissl_space/','nisslstack.vtk','atlas-seg','N',5);
% addpath(genpath('/Users/bhuo/Documents/GITHUB/Registration'))
% cd('/Users/bhuo/Dropbox (Mitra Lab)/Data and Analysis/Marmoset/Connectivity/M6328/Registered_Space_Outputs_v01/2D_nissl_space')
% filelist=filelsread('atlas-seg*F*.vtk','../',5);
% filenamepart='atlas-seg';
% modality='F';
% idloc=5;
cdir=pwd;
cd(vtkdir);
filelist=filelsread([filenamepart,'*',modality,'*.vtk'],'./',idloc);
[x,y,z,img1,title,names,spacing,origin] = read_vtk_image(filelist{1});
imgstack=zeros(size(img1,1),size(img1,2),length(filelist));
imgstack=cast(imgstack,'like',img1);
imgstack(:,:,1)=img1;
for f=2:length(filelist)
    disp(filelist{f})
    tic;[~,~,~,imgstack(:,:,f),~,~] = read_vtk_image(filelist{f},0);toc;
end
tic; vtkwrite(stackfilename,'structured_points',title,imgstack,'spacing',spacing(1),spacing(2),spacing(3),'origin',origin(1),origin(2),origin(3)); toc;
% tic; vtkwrite('test.vtk','structured_points',title,imgstack,'spacing',spacing(1),spacing(2),80); toc;
cd(cdir)