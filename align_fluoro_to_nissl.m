function align_fluoro_to_nissl(input_dir, nissl_pattern, fluoro_pattern, detailed_output_dir)
% for each floro image, this function computes a rigid alignment from the
% nearest nissl to it
% for reconstruction, this will need to be combined with transformation of
% the neighbor

addpath Functions/plotting

% get the target data
geometry_file = dir([input_dir '*.csv']);
fid = fopen([input_dir geometry_file(1).name],'rt');
line = fgetl(fid); % ignore the first line
csv_data = {};
is_nissl = [];
is_fluoro = [];
count = 0;
while 1
    line = fgetl(fid);
    if line == -1
        break
    end
    
    count = count + 1;
    
    % check if it matches the pattern
    if (regexp(line,regexptranslate('wildcard',nissl_pattern)))
        is_nissl(count) = 1;
        is_fluoro(count) = 0;
    elseif (regexp(line,regexptranslate('wildcard',fluoro_pattern)))
        is_fluoro(count) = 1;
        is_nissl(count) = 0;
    else
        disp([num2str(count)  ' ' line])
    end
    

    % process this line, splitting at commas
    csv_data(count,:) = strsplit(line,',');
    %     
end
fclose(fid);
files = csv_data(:,1);
nxJ0 = cellfun(@(x)str2num(x), csv_data(:,2:3));
dxJ0 = cellfun(@(x)str2num(x), csv_data(:,5:6));
zJ0 = cellfun(@(x)str2num(x), csv_data(:,10));

%%
% now pairwise registration
NtoF = repmat(eye(3),1,1,length(zJ0));
inds = zeros(1,length(zJ0));
for i = 1 : length(zJ0)
    if ~is_fluoro(i)
        continue
    end
    
    % find the nearest nissl
    cost = (zJ0(i) - zJ0).^2 + is_fluoro(:)*1e10;
    ind = find(  cost == min(cost) ,1,'first');
    inds(i) = ind;
    % load the images
    I = imread([input_dir files{ind}]);
    I = double(I)/255.0;
    dxI = dxJ0(ind,:);
    xI = (0:size(I,2)-1)*dxI(1);xI = xI - mean(xI);
    yI = (0:size(I,1)-1)*dxI(2);yI = yI - mean(yI);
    
    J = imread([input_dir files{i}]);
    J = double(J)/255.0;
    dxJ = dxJ0(i,:);
    xJ = (0:size(J,2)-1)*dxJ(1);xJ = xJ - mean(xJ);
    yJ = (0:size(J,1)-1)*dxJ(2);yJ = yJ - mean(yJ);
    
    
    
    A0 = eye(3);
    
    % we want an initial translation by center of mass
    % we will do this by thresholding on brightness and taking largest
    % connected component
    % maybe do EM with 3 compartments
    mask = getmask(I,1);

    % I is dark on light
    [XI,YI] = meshgrid(xI,yI);
    W = detect_padding_in_nissl(I);
    mask = mask & ~W;
    mask = double(mask) / sum(mask(:));
    comI = [sum(mask(:).*XI(:)),sum(mask(:).*YI(:))];
    if any(isnan(comI))
        comI = [0,0];
    end
    
    % J is light on dark
    [XJ,YJ] = meshgrid(xJ,yJ);
    [~,~,ext] = fileparts(fluoro_pattern);
    if strcmp(fluoro_pattern, ['*-F*' ext])
        mask = getmask(J,0);
    elseif strcmp(fluoro_pattern, ['*-IHC*' ext]) % dark on light
        mask = getmask(J,1);
    end
    mask =  double(mask)/sum(mask(:));
    comJ = [sum(mask(:).*XJ(:)),sum(mask(:).*YJ(:))];
    if any(isnan(comJ))
        comJ = [0,0];
    end
    A0(1:2,end) = comJ - comI;
    % my com calcs are not very robust
    % look at some images with artifacts and you'll see
    
    danfigure(51);
    imagesc(xI,yI,I);
    axis image
    hold on;
    scatter(comI(1),comI(2),'r','f')
    hold off;
    danfigure(52);
    imagesc(xJ,yJ,J);
    axis image
    hold on;
    scatter(comJ(1),comJ(2),'r','f')
    hold off;
    
    
    downs = [16,8,6,3]; % note that 3 takes us from 15 to 45, the res we're working at
    niter = 100 *ones(size(downs));
    niter(1) = 500;
    niter = [400,200,100,50];
    niter = [200,100,50,25];

    
%     eTfactor = 1e-2*5;
%     eLfactor = 1e-8*1;
%     NtoF(:,:,i) = slice_to_slice_rigid_alignment(xI,yI,I,xJ,yJ,J,A0,downs,niter,eTfactor,eLfactor);
    % this gauss newton version is much more numerically stable
    niter = 100;
    try
        NtoF(:,:,i) = slice_to_slice_rigid_alignment_GN_weight(xI,yI,I,xJ,yJ,J,A0,downs,niter);
    catch
        warning(['Could not calculate slice nissl to fluoro alignment for slice ' num2str(i) ', setting to identity']);
        NtoF(:,:,i) = eye(3);        
    end
end

%%
% now we have to save this data
save([detailed_output_dir 'NtoF.mat'],'NtoF','is_nissl','is_fluoro','files','inds')
