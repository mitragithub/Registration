% modified from Daniel Tward's transform.m and Adam Lin's transform_injection_atlas_full.m
% Inputs
function OUTPUTvol=transform_reg_to_atlas(datainfo,INPUTvol,vtkpath,atlas_vtk)
cd(vtkpath)
vtklist=filelsread(['registered_to_atlas_displacement_*',datainfo.staining,'*.vtk'],'~/',8);
[count2, b] = size(vtklist);
% initialize volume
[xA,yA,zA,atlas,title_,names,atlas_spacing] = read_vtk_image(atlas_vtk);
atlas_spacing=unique(atlas_spacing); % atlas should be isometric
C=size(INPUTvol,4);
OUTPUTvol = zeros([size(atlas)*round(atlas_spacing/datainfo.voxelsize),C]);

% transformed_path = strcat(output_dir, brain_no, '/', tracer, '/');

% imglist = dir(strcat(transformed_path, '*.tif'));
% %figure;
% [count, a] = size(imglist);
% disp(count);
% for i = 1 : count
%     imgpath = (strcat(imglist(i).folder, '/', imglist(i).name)); % image in registered space
%     imgNumber = split(imglist(i).name, '_');
%     imgNumber = imgNumber(4);
%     staining = split(imglist(i).name, '-');
%     staining = staining(2);
%%
for j = 1:count2
    Jhigh = squeeze(INPUTvol(:,:,j,:));
    if sum(sum(sum(Jhigh)))>0
        reconpath=vtklist{j}; %vtk list is sorted
        %         reconpath = (strcat(vtklist(j).folder, '/', vtklist(j).name)); % transformation vtk file
        %         reconpathNumber = split(vtklist(j).name, '_');
        %         reconpathNumber = erase(reconpathNumber(8),'.vtk');
        %         vtkstaining = split(vtklist(j).name, '-');
        %         vtkstaining = vtkstaining(2);
        %         if isequal(reconpathNumber,imgNumber) && isequal(vtkstaining,staining)
        %             disp(imgpath)
        disp(reconpath)
        %% transform.m
        
        % get its pixel size in um
        dxhigh = [datainfo.voxelsize,datainfo.voxelsize];
        
        
        % get the number of pixels
        nxhigh = [size(Jhigh,2),size(Jhigh,1)];
        % get the domain
        xhigh = (0:nxhigh(1)-1)*dxhigh(1);
        yhigh = (0:nxhigh(2)-1)*dxhigh(2);
        % make sure it is centered at 0
        xhigh = xhigh - mean(xhigh);
        yhigh = yhigh - mean(yhigh);
        % use meshgrid for sampling
        [XHIGH,YHIGH] = meshgrid(xhigh,yhigh);
        
        %% transform_injection_atlas_full.m
        % load transformation
        [x,y,z,I,title_,names] = read_vtk_image(reconpath);
        [X,Y,Z] = meshgrid(x,y,z);
        %     [X,Y] = meshgrid(x,y);
        
        AJphiJAxyX = I(:,:,1,1)+X;
        AJphiJAxyY = I(:,:,1,2)+Y;
        AJphiJAxyZ = I(:,:,1,3)+Z;
        %             image = imread(imgpath);
        
        %             [nx,ny] = size(image);
        %             test = zeros(160,228);
        
        %             F = griddedInterpolant({yA,xA,zA},atlas,'linear','nearest');
        %% transform.m
        % upsample the transformation
        F = griddedInterpolant({y,x},AJphiJAxyX,'linear','nearest');
        AJphiJXhigh = F(YHIGH,XHIGH);
        F = griddedInterpolant({y,x},AJphiJAxyY,'linear','nearest');
        AJphiJYhigh = F(YHIGH,XHIGH);
        F = griddedInterpolant({y,x},AJphiJAxyZ,'linear','nearest');
        AJphiJZhigh = F(YHIGH,XHIGH);
        %
        
        %     % apply transformation to your image
        %     transformed_Jhigh = zeros(size(Jhigh));
        %     for c = 1 : size(Jhigh,3)
        %         F = griddedInterpolant({yhigh,xhigh},Jhigh(:,:,c),'linear','none');
        %         transformed_Jhigh(:,:,c) = F(AJphiJYhigh,AJphiJXhigh);
        %     end
        %
        %
        %     transformed_Jhigh(isnan(transformed_Jhigh)) = 0;
        %     regatlasvol(:,:,j,:)=transformed_Jhigh;
        %                 transformed_Jhigh=im2uint8(transformed_Jhigh);
        %% transform_injection_atlas_full.m
        %             testdan = F(AJphiJAxyY,AJphiJAxyX,AJphiJAxyZ);
        %                 coordsX = AJphiJAxyX(image>0);
        %                 coordsY = AJphiJAxyY(image>0);
        %                 coordsZ = AJphiJAxyZ(image>0);
        for c = 1 : C
            image=Jhigh(:,:,c);
            if sum(sum(sum(image)))>0
                coordsX = AJphiJXhigh(image>0);
                coordsY = AJphiJYhigh(image>0);
                coordsZ = AJphiJZhigh(image>0);
                imagevals=image(image>0);
                
                for i = 1 : length(coordsX)
                    coordX = coordsX(i);
                    coordY = coordsY(i);
                    coordZ = coordsZ(i);
                    
                    indX = find( (xA - coordX).^2 == min((xA - coordX).^2),1,'first');
                    indY = find( (yA - coordY).^2 == min((yA - coordY).^2),1,'first');
                    indZ = find( (zA - coordZ).^2 == min((zA - coordZ).^2),1,'first');
                    %disp(indZ);
                    OUTPUTvol(indY,indX,indZ,c) = imagevals(i);
                end
            end
        end
    end
    
end