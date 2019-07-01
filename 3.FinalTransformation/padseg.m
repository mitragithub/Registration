function padseg(seg_path, output_dir)
    % LASTN = maxNumCompThreads(1)

    p = parpool('local',10);
    pctRunOnAll maxNumCompThreads(1)
    seglist = dir(strcat(seg_path, '*.mat'));
    parfor i = 1:length(seglist)
        segpath = (strcat(seglist(i).folder, '/', seglist(i).name));
        % load(segpath)
        % img = seg;
        img = padseg_aux(segpath);
        disp(segpath);
        [m n z] = size(img);
        p = 24000;
        q = 24000;
    
        img_pad = zeros([p q z]);
        for c = 1 : size(img,3)
            tmp = padarray(img(:,:,c), [floor((p - m)/2) floor((q - n)/2)], 0, 'post');
            img_pad(:,:,c) = padarray(tmp, [ceil((p - m)/2) ceil((q - n)/2)], 0, 'pre');
        end
    
        seg = img_pad;
        segresize_aux2(seg, strcat(output_dir, seglist(i).name));
        % save(strcat(output_dir, seglist(i).name),'seg','-v7.3')
    end
    delete(p)

end


% clear all;

% seglist = dir('/scratch/783UPDATE/reg_high_seg/*.mat');

% for i = 1:length(seglist)
%     segpath = (strcat(seglist(i).folder, '/', seglist(i).name));
%     load(segpath)
%     img = seg;
%     disp(segpath);
%     [m n z] = size(img);
%     p = 30000;
%     q = 30000;

%     img_pad = zeros([p q z]);
%     for c = 1 : size(img,3)
%         tmp = padarray(img(:,:,c), [floor((p - m)/2) floor((q - n)/2)], 0, 'post');
%         img_pad(:,:,c) = padarray(tmp, [ceil((p - m)/2) ceil((q - n)/2)], 0, 'pre');
%     end

%     seg = img_pad;
%     save(strcat('/scratch/783UPDATE/reg_high_seg_pad/', seglist(i).name),'seg','-v7.3')
% end