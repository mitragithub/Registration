function segresize(seg_path, img_path, output_dir)
    % LASTN = maxNumCompThreads(1)

    p = parpool('local',10);
    pctRunOnAll maxNumCompThreads(1)
    seglist = dir(strcat(seg_path, '*.mat'));
    imglist = dir(strcat(img_path, '*.tif'));
    [count, a] = size(imglist);

    parfor i = 1:count
        % LASTN = maxNumCompThreads(1)
        segpath = (strcat(seglist(i).folder, '/', seglist(i).name));
        imgpath = (strcat(imglist(i).folder, '/', imglist(i).name));
        disp(segpath);

        transformed_Jhigh = segresize_aux(segpath);
        seg = (transformed_Jhigh);
        img = imread(imgpath);
        [m n z] = size(img);
        seg = imresize(seg, [m n],'nearest');
        segresize_aux2(seg, strcat(output_dir, seglist(i).name))
        % save(strcat(output_dir, seglist(i).name), 'seg','-v7.3')
    end
    delete(p)

end



% seglist = dir('/scratch/783UPDATE/lowseg/*.mat'); %this is the registered low resolution segmentation
% imglist = dir('/scratch/783UPDATE/reg_high_tif/*.tif');


% for i = 1:length(seglist)
    
%         segpath = (strcat(seglist(i).folder, '/', seglist(i).name));
%         imgpath = (strcat(imglist(i).folder, '/', imglist(i).name));
%         disp(segpath);

%         load(segpath);%transformed_Jhigh, not a very good name, fix it.
%         seg = (transformed_Jhigh);
%         img = imread(imgpath);
%         [m n z] = size(img);
%         %disp(size(img))
%         seg = imresize(seg, [m n],'nearest');
%         save(strcat('/scratch/783UPDATE/reg_high_seg/', seglist(i).name), 'seg','-v7.3')
% end