function align_fluoro_to_nissl(input_dir, nissl_pattern, fluoro_pattern, detailed_output_dir)
% for each floro image, this function computes a rigid alignment from the
% nearest nissl to it
% for reconstruction, this will need to be combined with transformation of
% the neighbor

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
    % I is dark on light
    [XI,YI] = meshgrid(xI,yI);
    I_ = min(I,[],3);
    I_ = max(I_(:)) - I_;
    I_ = I_/sum(I_(:));
    comI = [sum(I_(:).*XI(:)),sum(I_(:).*YI(:))];
    if any(isnan(comI))
        comI = [0,0];
    end
    
    % J is light on dark
    if strcmp(fluoro_pattern, '*-F*.tif')
        [XJ,YJ] = meshgrid(xJ,yJ);
        J_ = max(J,[],3);
        J_ = J_ - min(J_(:));
    elseif strcmp(fluoro_pattern, '*-IHC*.tif') % dark on light
        [XJ,YJ] = meshgrid(xJ,yJ);
        J_ = min(J,[],3);
        J_ = max(J_(:)) - J_;
    end
    J_ = J_/sum(J_(:));
    comJ = [sum(J_(:).*XJ(:)),sum(J_(:).*YJ(:))];
    if any(isnan(comJ))
        comJ = [0,0];
    end
    A0(1:2,end) = comJ - comI;
    % my com calcs are not very robust
    % look at some images with artifacts and you'll see
    
    danfigure(1);
    imagesc(xI,yI,I);
    axis image
    hold on;
    scatter(comI(1),comI(2),'r','f')
    hold off;
    danfigure(2);
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
        NtoF(:,:,i) = slice_to_slice_rigid_alignment_GN(xI,yI,I,xJ,yJ,J,A0,downs,niter);
    catch
        warning(['Could not calculate slice nissl to fluoro alignment for slice ' num2str(i)]);
        keyboard
    end
end

%%
% now we have to save this data
save([detailed_output_dir 'NtoF.mat'],'NtoF','is_nissl','is_fluoro','files','inds')
