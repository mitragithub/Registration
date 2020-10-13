addpath(genpath('/Users/bhuo/Documents/GITHUB/Registration'))
cd('/Users/bhuo/Dropbox (Mitra Lab)/Data and Analysis/Marmoset/Connectivity/M6328/Registered_Space_Outputs_v01/2D_nissl_space')
filelist=filelsread('atlas-seg*F*.vtk','../',5);
 [x,y,z,img1,title,names] = read_vtk_image(filelist{1});
 imgstack=zeros(size(img1,1),size(img1,2),length(filelist));
 imgstack=cast(imgstack,'like',img1);
 imgstack(:,:,1)=img1;
for f=2:length(filelist)
    disp(filelist{f})
    tic;[x,y,z,imgstack(:,:,f),title_,names] = read_vtk_image(filelist{f},0);toc;
end
saveimgstack(imgstack,'