function padtif(img_path, output_dir, output_jp2)
    % LASTN = maxNumCompThreads(1)

    parpool('local',10);
    pctRunOnAll maxNumCompThreads(1)
    imglist = dir(strcat(img_path,'/*.tif'));
    parfor i = 1:length(imglist)
        % LASTN = maxNumCompThreads(1)

        imgpath = (strcat(imglist(i).folder, '/', imglist(i).name));
        img = imread(imgpath);
        disp(imgpath);
        [m n z] = size(img);
        p = 24000;
        q = 24000;
    
        img_pad = uint8(zeros([p q z]));
        for c = 1 : size(img,3)
            tmp = padarray(img(:,:,c), [floor((p - m)/2) floor((q - n)/2)], 255, 'post');
            img_pad(:,:,c) = padarray(tmp, [ceil((p - m)/2) ceil((q - n)/2)], 255, 'pre');
        end
        imwrite(img_pad, strcat(output_dir, imglist(i).name),'Compression', 'none')
        % disp(['sh kducomp_single.sh "' strcat(output_dir, imglist(i).name) '" ' output_jp2])
        system(['sh kducomp_single.sh "' strcat(output_dir, imglist(i).name) '" ' output_jp2]);
    end
end
% clear all;

% imglist = dir('/scratch/783UPDATE/reg_high_tif/*.tif');

% parfor i = 1:length(imglist)
%     imgpath = (strcat(imglist(i).folder, '/', imglist(i).name));
%     img = imread(imgpath);
%     disp(imgpath);
%     [m n z] = size(img);
%     p = 30000;
%     q = 30000;

%     img_pad = uint8(zeros([p q z]));
%     for c = 1 : size(img,3)
%         tmp = padarray(img(:,:,c), [floor((p - m)/2) floor((q - n)/2)], 255, 'post');
%         img_pad(:,:,c) = padarray(tmp, [ceil((p - m)/2) ceil((q - n)/2)], 255, 'pre');
%     end
%     imwrite(img_pad, strcat('/scratch/783UPDATE/reg_high_tif_pad/', imglist(i).name),'Compression', 'none')
% end