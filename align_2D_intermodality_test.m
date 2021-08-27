function align_2D_intermodality_test(input_dir, pattern0, pattern1, detailed_output_dir, OPT)

% this test is for marmoset data
% in test I'm trying to improve the myelin maps by matching based on masks
% 1. putting somothing into gmm was helpful
% 2. switching to mutual information was helpful
% 3. remaining problem is bad mask generation in posterior, I added an
% extra compartment to masking, seems okay
%

% rigidly align 2D intermodality image slices using polynomial contrast 
% prediction and missing data estimation
% For each image matching pattern1, the closest image matching pattern0 is
% identified.  Rigid registration is then applied to align the pattern0
% image to the pattern1 image, and each registaration is saved as an affine
% matrix.  Registration is initialized with a translation by identifying an 
% image mask and finding its center of mass.
%
% arguments:
% input_dir:  directory containing 2D low resolution images and geometry 
%             csv information files
% pattern0:   glob pattern (e.g. with wildcards) to identify files to match
%             from (i.e. moving image).  Typically this is a Nissl image.
% pattern1:   glob pattern (e.g. with wildcards) to identify files to match
%             to (i.e. fixed image).  Typically this is a fluorescence 
%             image.
% detailed_output_dir: location where outputs will be saved.
% OPT:        Other optional inputs as a structure. specific usage
%             described below. Experimental
% outputs:
%             returns nothing but writes out a .mat file saving rigid
%             transformation matrices in the detailed_output_dir
%

addpath Functions/plotting

% keyboard
skip = 1; % for debugging set to a big number
if nargin < 5
    OPT = struct;
end
if isfield(OPT,'method')
    method = lower(OPT.method);
else
    method = 'mi';
end
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
%         disp([num2str(count)  ' ' line]) % e.g. if I have nissl myelin
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
mapcount = 0;
for i = 1 : length(zJ0)
    if ~is_1(i)
        continue
    end
    mapcount = mapcount + 1;
    if mod(mapcount,skip)
        continue
    end
    
    
    
    % find the nearest pattern 0
    cost = (zJ0(i) - zJ0).^2 - is_0(:)*1e10;
    ind = find(  cost == min(cost) ,1,'first');
    inds(i) = ind;
    % load the images
    I = imread([input_dir files{ind}]);
    if isa(I,'uint8')
        I = double(I)/255.0;    
    elseif isa(I,'uint16')
        I = double(I);
        for c = 1 : size(I,3)
            I(:,:,c) = tiedrank(I(:,:,c));
        end
        I = I / max(I(:));
    end
    dxI = dxJ0(ind,:);
    xI = (0:size(I,2)-1)*dxI(1);xI = xI - mean(xI);
    yI = (0:size(I,1)-1)*dxI(2);yI = yI - mean(yI);
    
    J = imread([input_dir files{i}]);
    if isa(J,'uint8')
        J = double(J)/255.0;    
    elseif isa(J,'uint16')
        
%         J = double(J);
%         for c = 1 : size(I,3)
%             J(:,:,c) = tiedrank(J(:,:,c));
%         end
%         J = J / max(J(:));
        
        J = double(J);
        J_ = J(:);
        J_ = tiedrank(J_);
        J = reshape(J_,size(J));
        J = J / max(J(:));
        
    end
    
    dxJ = dxJ0(i,:);
    xJ = (0:size(J,2)-1)*dxJ(1);xJ = xJ - mean(xJ);
    yJ = (0:size(J,1)-1)*dxJ(2);yJ = yJ - mean(yJ);
    
    
    A0 = eye(3);
    
    % we want an initial translation by center of mass
    % we will do this by thresholding on brightness and taking largest
    % connected component
    

    
    
    
    % I
    % separate it into two components, we want mask to indicate foreground
%     mask = kmeans(reshape(I,[size(I,1)*size(I,2),size(I,3)]),2);
%     mask = reshape(mask(:,:,1),size(I,1),size(I,2));
    nblur = 3;
    e = 0.5;

    ngmm = 50;
    % two components, light and dark
    mu = [[1,1,1];[0.5,0.5,0.5];[0,0,0];]';
    ncomp = size(mu,2);
    mu = reshape(mu,[1,3,ncomp]);
    p = ones(1,ncomp)/ncomp;
    Sigma = repmat(eye(3),[1,1,ncomp]);
    X = reshape(I,[size(I,1)*size(I,2),size(I,3)]);
    for it = 1 : ngmm
        % given mu and sigma, calculate probs
        % subtract mean
        X0 = X-mu;
        Sigmai = zeros(3,3,ncomp);
        detSigma = zeros(1, ncomp);
        for c = 1 : ncomp
            Sigmai(:,:,c) = inv(Sigma(:,:,c));
            detSigma(c) = det(Sigma(:,:,c));
        end
        % quadratic form
        Q = zeros(size(X,1),size(Sigma,3));
        for c = 1 : size(Sigma,3)
            Q(:,c) = sum( (X0(:,:,c) * Sigmai(:,:,c)) .* X0(:,:,c), 2);
        end
        % likelihood
        l = exp(-Q/2.0)./sqrt(2.0*pi*detSigma);
        l = (l.*p)./(sum(l.*p,2) + max(sum(l.*p,2))*1e-12);
       
        n = sum(l,1);
        p =  n/sum(n);
        
        % blur l
        L = reshape(l,size(I,1),size(I,2),ncomp);
        for a = 1 : nblur
        L = (1-e).*L + e/4.*(circshift(L,[1,0]) + circshift(L,[-1,0]) + circshift(L,[0,1]) + circshift(L,[0,-1]));
        end
        l = reshape(L,[],ncomp);

        w = l./n;

        
        % given probs calculate mu and sigma
        mu = sum(reshape(X,size(X,1),size(X,2),1).*reshape(w,size(w,1),1,size(w,2)),1);
        X0 = X - mu;
        for c = 1 : size(Sigma,3)
            Sigma(:,:,c) = X0(:,:,c)' * (X0(:,:,c).*w(:,c));
        end
        
        % blur mu and sigma
        if ~mod(it,10) || it == ngmm
        danfigure(55);
        imagesc(reshape(l(:,2),size(I,1),size(I,2)))
        axis image
        drawnow
        end
        
    end
    mask = reshape(l(:,1),size(I,1),size(I,2)); % brightest for nissl is background
    
    
    
    % look at the edges
    if mean( [mask(1,:),mask(end,:),mask(:,1)',mask(:,end)'] ) > 0.5
        mask = 1 - mask;
    end
    maskI = mask;
    [XI,YI] = meshgrid(xI,yI);
    mask = double(mask) / sum(mask(:));
    comI = [sum(mask(:).*XI(:)),sum(mask(:).*YI(:))];
    if any(isnan(comI))
        comI = [0,0];
    end
    
    % J 
    % two components, light and dark
    mu = [[1,1,1];[0,0,0]]';
    ncomp = size(mu,2);
    mu = reshape(mu,[1,3,ncomp]);
    p = ones(1,ncomp)/ncomp;
    Sigma = repmat(eye(3),[1,1,ncomp]);
    X = reshape(J,[size(J,1)*size(J,2),size(J,3)]);
    for it = 1 : ngmm
        % given mu and sigma, calculate probs
        % subtract mean
        X0 = X-mu;
        Sigmai = zeros(3,3,ncomp);
        detSigma = zeros(1, ncomp);
        for c = 1 : ncomp
            Sigmai(:,:,c) = inv(Sigma(:,:,c));
            detSigma(c) = det(Sigma(:,:,c));
        end
        % quadratic form
        Q = zeros(size(X,1),size(Sigma,3));
        for c = 1 : size(Sigma,3)
            Q(:,c) = sum( (X0(:,:,c) * Sigmai(:,:,c)) .* X0(:,:,c), 2);
        end
        % likelihood
        l = exp(-Q/2.0)./sqrt(2.0*pi*detSigma);
        l = (l.*p)./(sum(l.*p,2) + max(sum(l.*p,2))*1e-12);
       
        n = sum(l,1);
        p =  n/sum(n);
        
        % blur l
        L = reshape(l,size(J,1),size(J,2),2);
        for a = 1 : nblur
        L = (1-e).*L + e/4.*(circshift(L,[1,0]) + circshift(L,[-1,0]) + circshift(L,[0,1]) + circshift(L,[0,-1]));
        end
        l = reshape(L,[],2);
        
        w = l./n;

        
        % given probs calculate mu and sigma
        mu = sum(reshape(X,size(X,1),size(X,2),1).*reshape(w,size(w,1),1,size(w,2)),1);
        X0 = X - mu;
        for c = 1 : size(Sigma,3)
            Sigma(:,:,c) = X0(:,:,c)' * (X0(:,:,c).*w(:,c));
        end
        if ~mod(it,10) || it == ngmm
        danfigure(56);
        imagesc(reshape(l(:,2),size(J,1),size(J,2)))
        axis image
        drawnow
        end
        
    end
    mask = reshape(l(:,2),size(J,1),size(J,2)); % bright

    
    % look at the edges
    if mean( [mask(1,:),mask(end,:),mask(:,1)',mask(:,end)'] ) > 0.5
        mask = 1 - mask;
    end
    maskJ = mask;
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
    downs = [4,2,1];
%     downs = [2,1];
%     downs = [1];
    
    % this gauss newton version is much more numerically stable
    disp(i)
    niter = 100;
%     niter = 300;
%     niter = [200,100,50];
    niter = [300,150,75];
    niter = [1000,100,50];
    niter = [1000,200,100];
    
    
    try
        
        if strcmp(method,'mi')
            bins = linspace(0,1,10);
            sigma = 0.2;
            e = 1e4; % before normalizing
    %         e = 2e5;

            sigma = 0.1;
            e = 5e2;
            % pca?        

            if strcmp(pattern1,'*-F*')
                J(:,:,2) = 0.0;
            end


            zero_to_one(:,:,i) = slice_to_slice_rigid_alignment_MI(xI,yI,I,xJ,yJ,J,A0,downs,niter,e,bins,sigma);
        elseif strcmp(method,'em')
            zero_to_one(:,:,i) = slice_to_slice_rigid_alignment_GN_weight(xI,yI,I,xJ,yJ,J,A0,downs,niter);
        elseif strcmp(method,'test1')
            % note this just uses masks
            OPT.nocontrast = 1;
            OPT.noweight = 0;
            OPT.sigmaM = 0.25;
            OPT.sigmaA = OPT.sigmaM*3;
            zero_to_one(:,:,i) = slice_to_slice_rigid_alignment_GN_weight(xI,yI,maskI,xJ,yJ,maskJ,A0,downs,niter,0.2,OPT);
        elseif strcmp(method,'test2')
            
            
            % we are going to replace I with some features
            % first pca
            X = reshape(I,[],size(I,3));
            Xbar = mean(X,1);
            X0 = X - Xbar;
            Sigma = X0'*X0/size(X0,1);
            [v,d] = eig(Sigma,eye(size(I,3)),'chol');
            % these should be sorted increasing
            I_ = reshape(X*v(:,end),size(I,1),size(I,2));
            I__ = reshape(X*v(:,end-1),size(I,1),size(I,2));
            [I_x,I_y] = gradient(I_);
            nGI = sqrt(I_x.^2 + I_y.^2);
            I_ = normcdf((I_-mean(I_(:)))/std(I_(:)));
            I__ = normcdf((I__-mean(I__(:)))/std(I__(:)));
            nGI = normcdf((nGI - mean(nGI(:)))/std(nGI(:)));
            
            
            zero_to_one(:,:,i) = slice_to_slice_rigid_alignment_GN_weight(xI,yI,cat(3,I_,I__,nGI),xJ,yJ,J,A0,downs,niter,0.2,OPT);

        end
        
    catch
        
        warning(['Could not calculate slice nissl to fluoro alignment for slice ' num2str(i) ', setting to identity']);
        zero_to_one(:,:,i) = eye(3);        
    end
end

%%
% now we have to save this data
save([detailed_output_dir 'zero_to_one_' strrep(pattern0,'*','x') '_to_' strrep(pattern1,'*','x') '.mat'],'zero_to_one','is_0','is_1','files','inds')
