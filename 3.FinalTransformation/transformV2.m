function transformV2(geometry_path, img_path, recon_path, output_dir)
    res = 0.46 * 2;
    % parpool('local',16);

    geometry_file = dir([geometry_path '*lossless.csv']);
    % disp([geometry_path '*.csv'])
    fid = fopen([geometry_path geometry_file(1).name],'rt');
    line = fgetl(fid); % ignore the first line

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

    for f = 1 : length(files)
        [directory,fname,ext] = fileparts(files{f});
        lossless_file = dir([img_path fname '_lossless.jp2']);
        % disp([lossless_file(1).folder '/' lossless_file(1).name])
        % % disp([img_path fname '_lossless.jp2'])
        J = imread([lossless_file(1).folder '/' lossless_file(1).name]);
        disp(size(J))

        if strcmp(class(J),'uint8')
            J = double(J)/255.0;
        end
        J = double(J);
        J = imresize(J,0.5,'nearest');

         %%
        % transform this image to registered space (2d transforms, identity in z)
        % this involves two steps
        % 1. upsample the transform to resolution res
        % 2. resample the data
        % 
        % load the transform
        try
            [xTJ,yTJ,zTJ,DeltaTJ,title_,names] = read_vtk_image([recon_path 'registered_to_input_displacement_' fname '.vtk']);
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
        disp([output_dir fname '.tif'])
        imwrite(JT, [output_dir fname '.tif'])

        break
    end




    % imglist = dir(strcat(img_path, '*_lossless.jp2'));
    % vtklist = dir(strcat(recon_path, 'input_to_registered_displacement*.vtk'));
    % [count, a] = size(imglist);
    % % disp(img_path)
    % % disp(recon_path)
    % % disp(size(imglist))
    % % disp(size(vtklist))
    % for i = 1 : 100 %count
    %     reconpath = (strcat(vtklist(i).folder, '/', vtklist(i).name));
    %     imgpath = (strcat(imglist(i).folder, '/', imglist(i).name));
    %     disp(imgpath)
    %     disp(reconpath)
    %     Jhigh = imread(imgpath);
    %     Jhigh = im2double(Jhigh);


    %     % !FIXME due to space and memory limit, use lossy resolution(dowmsaple 2 times)
    %     Jhigh = imresize(Jhigh,0.5,'nearest');
    %     % get its pixel size in um
    %     dxhigh = [0.46 * 2, 0.46 * 2]; % just guessing here, use your actual pixel size

        
    %     % get the number of pixels
    %     nxhigh = [size(Jhigh,2),size(Jhigh,1)];
    %     % get the domain
    %     xhigh = (0:nxhigh(1)-1)*dxhigh(1);
    %     yhigh = (0:nxhigh(2)-1)*dxhigh(2);
    %     % make sure it is centered at 0
    %     xhigh = xhigh - mean(xhigh);
    %     yhigh = yhigh - mean(yhigh);
    %     % use meshgrid for sampling
    %     [XHIGH,YHIGH] = meshgrid(xhigh,yhigh);


    %     % load the transformation
    %     % [x, y, z, AJphiJAxyX, AJphiJAxyY] = transform_aux(reconpath);
    %     [x,y,z,I,title,names] = read_vtk_image(reconpath);
    %     [X,Y] = meshgrid(x,y);

    %     AJphiJAxyX = I(:,:,1,1)+X;
    %     AJphiJAxyY = I(:,:,1,2)+Y;
    %     % upsample the transformation
    %     F = griddedInterpolant({y,x},AJphiJAxyX,'linear','nearest'); % note x,y are 
    %     %loaded from the mat file, note always use order y,x
    %     AJphiJXhigh = F(YHIGH,XHIGH);
    %     F = griddedInterpolant({y,x},AJphiJAxyY,'linear','nearest');
    %     AJphiJYhigh = F(YHIGH,XHIGH);

    %     % apply transformation to your image
    %     transformed_Jhigh = zeros(size(Jhigh));
    %     for c = 1 : size(Jhigh,3)
    %         F = griddedInterpolant({yhigh,xhigh},Jhigh(:,:,c),'linear','none');
    %         transformed_Jhigh(:,:,c) = F(AJphiJYhigh,AJphiJXhigh);
    %     end

    %     transformed_Jhigh(isnan(transformed_Jhigh)) = 255;
    %     transformed_Jhigh=im2uint8(transformed_Jhigh);
    %     imwrite(transformed_Jhigh, strcat(output_dir, imglist(i).name(1:end-13),'.tif'))


    % end
end


% for i = 1:50
%     %try
%         reconpath = (strcat(reconlist(i).folder, '/', reconlist(i).name));
%         imgpath = (strcat(imglist(i).folder, '/', imglist(i).name));
%         disp(imgpath)
%         %if exist(strcat('/scratch/783UPDATE/reg_high_tif/', imglist(i).name(1:end-13),'.tif'),'file')==2
%            % disp('we got it!')
%             %continue
%         %end
%         % load your high res image for the "f"th slice
%         Jhigh = imread(imgpath);
%         Jhigh = im2double(Jhigh);

%         Jhigh = imresize(Jhigh,0.5,'nearest');
%         % get its pixel size in um
%         dxhigh = [0.46 * 2, 0.46 * 2]; % just guessing here, use your actual pixel size
%         % get the number of pixels
%         nxhigh = [size(Jhigh,2),size(Jhigh,1)];
%         % get the domain
%         xhigh = (0:nxhigh(1)-1)*dxhigh(1);
%         yhigh = (0:nxhigh(2)-1)*dxhigh(2);
%         % make sure it is centered at 0
%         xhigh = xhigh - mean(xhigh);
%         yhigh = yhigh - mean(yhigh);
%         % use meshgrid for sampling
%         [XHIGH,YHIGH] = meshgrid(xhigh,yhigh);

%         % load the transformation
%         load(reconpath)
%         % upsample the transformation
%         F = griddedInterpolant({y,x},AJphiJAxyX,'linear','nearest'); % note x,y are 
%         %loaded from the mat file, note always use order y,x
%         AJphiJXhigh = F(YHIGH,XHIGH);
%         F = griddedInterpolant({y,x},AJphiJAxyY,'linear','nearest');
%         AJphiJYhigh = F(YHIGH,XHIGH);

%         % apply transformation to your image
%         transformed_Jhigh = zeros(size(Jhigh));
%             for c = 1 : size(Jhigh,3)
%                 F = griddedInterpolant({yhigh,xhigh},Jhigh(:,:,c),'linear','none');
%                 transformed_Jhigh(:,:,c) = F(AJphiJYhigh,AJphiJXhigh);
%             end

%         transformed_Jhigh(isnan(transformed_Jhigh)) = 255;
%         transformed_Jhigh=im2uint8(transformed_Jhigh);
%         imwrite(transformed_Jhigh, strcat('/scratch/783UPDATE/reg_high_tif/', imglist(i).name(1:end-13),'.tif'))
%     %catch
%         %disp(strcat(imgpath,' not working!'))
%     %end
% end

