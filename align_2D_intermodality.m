function align_2D_intermodality(input_dir, pattern0, pattern1, detailed_output_dir)
% align 2D intermodality images
% map images matching pattern 0 to the closest image matching pattern 1

addpath Functions/plotting

% get the target data
geometry_file = dir([input_dir '*.csv']);
fid = fopen([input_dir geometry_file(1).name],'rt');
line = fgetl(fid); % ignore the first line
csv_data = {};
is_0 = [];
is_1 = [];
count = 0;
while 1
    line = fgetl(fid);
    if line == -1
        break
    end
    
    count = count + 1;
    
    % check if it matches the pattern
    if (regexp(line,regexptranslate('wildcard',pattern0)))
        is_0(count) = 1;
        is_1(count) = 0;
        is_other(count) = 0;
    elseif (regexp(line,regexptranslate('wildcard',pattern1)))
        is_1(count) = 1;
        is_0(count) = 0;
        is_other(count) = 0;
    else
        disp([num2str(count)  ' ' line]) % e.g. if I have nissl myelin
        is_other(count) = 1;
        is_0(count) = 0;
        is_1(count) = 0;
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
% now pairwise registration, loop through all the pattern 1
zero_to_one = repmat(eye(3),1,1,length(zJ0));
inds = zeros(1,length(zJ0));
for i = 1 : length(zJ0)
    if ~is_1(i)
        continue
    end
    
    
    % find the nearest pattern 0
    cost = (zJ0(i) - zJ0).^2 - is_0(:)*1e10;
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
    
    % start with a rank transform
    for c = 1 : size(I,3)
        I(:,:,c) = tiedrank(I(:,:,c));
    end
    for c = 1 : size(J,3)
        J(:,:,c) = tiedrank(J(:,:,c));
    end
    I = I/max(I(:));
    J = J/max(J(:));
    
    
    % I
    % separate it into two components, we want mask to indicate foreground
    mask = kmeans(reshape(I,[size(I,1)*size(I,2),size(I,3)]),2);
    mask = reshape(mask(:,:,1),size(I,1),size(I,2));
    % look at the edges
    if mean( [mask(1,:),mask(end,:),mask(:,1)',mask(:,end)'] ) > 0.5
        mask = 1 - mask;
    end
    [XI,YI] = meshgrid(xI,yI);
    mask = double(mask) / sum(mask(:));
    comI = [sum(mask(:).*XI(:)),sum(mask(:).*YI(:))];
    if any(isnan(comI))
        comI = [0,0];
    end
    
    % J 
    mask = kmeans(reshape(J,[size(J,1)*size(J,2),size(J,3)]),2);
    mask = reshape(mask(:,:,1),size(J,1),size(J,2));
    % look at the edges
    if mean( [mask(1,:),mask(end,:),mask(:,1)',mask(:,end)'] ) > 0.5
        mask = 1 - mask;
    end
    [XJ,YJ] = meshgrid(xJ,yJ);
    mask =  double(mask)/sum(mask(:));
    comJ = [sum(mask(:).*XJ(:)),sum(mask(:).*YJ(:))];
    if any(isnan(comJ))
        comJ = [0,0];
    end
    A0(1:2,end) = comJ - comI;

    
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
    niter = [200,100,50,25];
    downs = [8,6,3];
    niter = [100,50,25];

    
    % this gauss newton version is much more numerically stable
    niter = 100;
    
    try
        zero_to_one(:,:,i) = slice_to_slice_rigid_alignment_GN_weight(xI,yI,I,xJ,yJ,J,A0,downs,niter);
    catch
        warning(['Could not calculate slice nissl to fluoro alignment for slice ' num2str(i) ', setting to identity']);
        zero_to_one(:,:,i) = eye(3);        
    end
end

%%
% now we have to save this data
save([detailed_output_dir 'zero_to_one_' strrep(pattern0,'*','x') '_to_' strrep(pattern1,'*','x') '.mat'],'zero_to_one','is_0','is_1','files','inds')
