function atlas_free_rigid_alignment(target_dir, pattern, output_dir, r, downs, niter, e, skip_thick, load_initializer)

% atlas free slice alignment 
% rigid alignment of a stack
% each slice is iteratively aligned to a weighted average of its neighbors
% this can be used as a standalone program
% or as an initial guess atlas based slice alignment
%
% TODO, deal with slice spacing in physical units
% TODO, deal with slices left out
% skip_thick is a number that we don't update if it is bigger than this
% going to make a large change here and switch to GN optimization

% keyboard

addpath Functions/plotting
addpath Functions/downsample
addpath Functions/gradient
addpath Functions/frame2Gif
addpath Functions/vtk

if nargin < 4
    r = 25;
end
if nargin < 5
    downs = [32,16,8,4];
end
if nargin < 6
    niter = 200;
end
if nargin < 7
    % step sizes
    e = 0.5;

end

if nargin < 8
    skip_thick = -1;
end
if nargin < 9
    load_initializer=0;
end







geometry_file = dir([target_dir '*.csv']);
fid = fopen([target_dir geometry_file(1).name],'rt');
line = fgetl(fid); % ignore the first line
% it should say
% filename, nx, ny, nz, dx, dy, dz, x0, y0, z0

csv_data = {};
count = 0;
while 1
    line = fgetl(fid);
    if line == -1
        break
    end
    % check if it matches the pattern
    if isempty(regexp(line,regexptranslate('wildcard',pattern)))
        continue
    end
    count = count + 1;
    % process this line, splitting at commas
    csv_data(count,:) = strsplit(line,',');
    %     
end
fclose(fid);
files = csv_data(:,1);
nxJ0 = cellfun(@(x)str2num(x), csv_data(:,2:3));
dxJ0 = cellfun(@(x)str2num(x), csv_data(:,5:6));
zJ0 = cellfun(@(x)str2num(x), csv_data(:,10));
dz0 = cellfun(@(x)str2num(x), csv_data(:,7));



% keyboard



n = length(files);


% reconstruction params
theta = zeros(1,n) ;
tx = zeros(1,n);
ty = zeros(1,n);

% AJ is the transform from RECONSTRUCTED to OBSERVED
% images deform with inverse, so we will use this to transform OBSERVED to
% RECONSTRUCTED
AJ = zeros(3,3,length(files));
% instead of zero let's initialize it
if exist([output_dir 'initializer_A.mat'],'file') && load_initializer
    vars = load([output_dir 'initializer_A.mat']);
    AJ = vars.AJ;

    % get angle from x,y with correct sign
    theta = squeeze(atan2(AJ(1,2,:),AJ(1,1,:)))';
    
    tx = squeeze(AJ(1,3,:))';
    ty = squeeze(AJ(2,3,:))';   
else
    for i = 1 : size(AJ,3)
        AJ(:,:,i) = eye(3);
    end
end

Cost = zeros(1,10000);

itercount = 0;

if length(niter) == 1
    niter = repmat(niter,length(downs));
end

% plotting
qlim = [0.001,0.999];
nplot = 5;
for downcount = 1 : length(downs)
    down = downs(downcount);

    
    % loop through the files
    nxJ = zeros(length(files),2);
    for i = 1 : n
        
        if i == round(n/2) && contains(target_dir,'human')
%             keyboard
        end
        J_ = imread([target_dir files{i}]);
        J_ = double(J_)/255.0;
        J_ = mean(J_,3);
        xJ{i} = (1 : size(J_,2))*dxJ0(i,1); 
        xJ{i} = xJ{i} - mean(xJ{i});
        yJ{i} = (1 : size(J_,1))*dxJ0(i,2); 
        yJ{i} = yJ{i} - mean(yJ{i});        
        
        % padding values
        val(i) = mode(J_(J_(:)~=1));
        if isnan(val(i))
            val(i) = val(i-1); % entirely missing slice?
        end
        J_(J_==1) = val(i); % take care of white background (padding that may have been inserted in other steps)
        
        
        % for thick images, just ignore them
        if dz0(i) > skip_thick && skip_thick > 0
            if i > 1
                J_ = ones(size(J_))*val(i-1);
            else
                J_ = zeros(size(J_));
            end
        end
        
        % take care of very dark bits
        thresh = 0.01;
        if contains(target_dir,'human')
            thresh = 0.15; % 0.15 seems to get the dark border for human data
        end
        W_ = J_>=thresh;
        J_(J_<thresh) = val(i);
        
        
        % downsample
        [~,~,W{i}] = downsample2D(xJ{i},yJ{i},W_,down*[1 1]);
        [xJ{i},yJ{i},J{i}] = downsample2D(xJ{i},yJ{i},J_.*W_,down*[1 1]);

        J{i} = J{i} ./ W{i};
        J{i}(isnan(J{i})) = val(i);
        
        % sometimes images will just be one pixel and downsampling doesn't work
        if size(J{i},1) < 2
            J{i} = [J{i};J{i}];
            yJ{i} = [yJ{i}(1),yJ{i}(1)+1];
        end
        if size(J{i},2) < 2
            J{i} = [J{i},J{i}];
            xJ{i} = [xJ{i}(1),xJ{i}(1)+1];
        end
        % no need to standardize
%         J{i} = ( J{i} - mean(J{i}(:)) ) / std(J{i}(:));
        
        nxJ(i,:) = [size(J{i},2),size(J{i},1)];
        dxJ(i,:) = [xJ{i}(2)-xJ{i}(1), yJ{i}(2)-yJ{i}(1)];
%         dxJ(i,:) = dxJ0(i,:)*down;
        
        % build some weights to avoid boundary effects
        try
            n_boundary = round(min([size(J{i},1),size(J{i},2)])*0.25); % about 20 percent boundary  
            tmp = build_boundary_weight([size(J{i},1),size(J{i},2)],n_boundary);
            W{i} = W{i}.*tmp;
        catch
            W{i} = W{i}.*ones([size(J{i},1),size(J{i},2)]);
        end
        
        % maybe contrast normalization will be useful
        if contains(target_dir,'human')
%             Jmean = sum(J{i}(:).*W{i}(:)) / sum(W{i}(:));
%             J2mean = sum(J{i}(:).^2.*W{i}(:)) / sum(W{i}(:));
%             
%             Jstd = J2mean - Jmean.^2;
%             Jstd(Jstd < 0) = 0;
%             Jstd = sqrt(Jstd);
%             Jz = (J{i} - Jmean)/Jstd;
%             J{i} = normcdf(Jz);
%             val(i) = normcdf((val(i) - Jmean)/Jstd);
        end
        
    end
    nxI = [max(nxJ,[],1),length(files)];
    % last component doesn't matter, display as a square
    dxI = [dxJ(1,:),1/n*nxI(1)*dxJ0(1,1)*down]; % common to all slices
    I0 = zeros(nxI(2),nxI(1),nxI(3));
    W0 = zeros(nxI(2),nxI(1),nxI(3));
    

    
    xI = (0 : nxI(1)-1)*dxI(1);xI = xI - mean(xI);
    yI = (0 : nxI(2)-1)*dxI(2);yI = yI - mean(yI);
    zI = (0 : nxI(3)-1)*dxI(3);zI = zI - mean(zI);
    [XI_,YI_] = meshgrid(xI,yI);
    

    
    
    %%
    % weights for averaging
    tmp = exp(-[-r:r].^2/2/(r/3)^2);
    tmp = tmp / sum(tmp);
    CC = tmp;
%     if ~contains(target_dir,'human')
    CC(([-r:r]==0)) = 0;
%     end
    CC = CC/sum(CC);
    

    
    
    %%
    for it = 1 : niter(downcount)        
        
        % load data again into I
        for i = 1 : n
            A = AJ(:,:,i);
            Xs = A(1,1)*XI_ + A(1,2)*YI_ + A(1,3);
            Ys = A(2,1)*XI_ + A(2,2)*YI_ + A(2,3);
            F = griddedInterpolant({yJ{i},xJ{i}},J{i},'linear','none');
            tmp = F(Ys,Xs);
            tmp(isnan(tmp)) = val(i);
            I0(:,:,i) = tmp;
            
            F = griddedInterpolant({yJ{i},xJ{i}},W{i},'linear','none');
            tmp = F(Ys,Xs);
            tmp(isnan(tmp)) = 0;
            W0(:,:,i) = tmp;
        end

    
    
        % first pad the image, and get a weighted average
        npad = r;
        I0p = padarray(I0.*W0,[0,0,npad],'both','symmetric');
        I = convn(I0p,reshape(CC,1,1,[]),'same');
        I = I(:,:,npad+1:end-npad);
        W0p = padarray(W0,[0,0,npad],'both','symmetric');
        W_ = convn(W0p,reshape(CC,1,1,[]),'same');
        W_ = W_(:,:,npad+1:end-npad);
        I = I ./ W_;
        
        
        % compute the image
        for i = 1 : n % loop over slices
            for sliceit = 1 : 1 % extra iterations
                
            if dz0(i) > skip_thick && skip_thick>0
                continue
            end
            % deform each slice
            Ai = inv(AJ(:,:,i));
            A = AJ(:,:,i);
            [XJ,YJ] = meshgrid(xJ{i},yJ{i});
            Xs = Ai(1,1)*XJ + Ai(1,2)*YJ  + Ai(1,3);
            Ys = Ai(2,1)*XJ + Ai(2,2)*YJ + Ai(2,3);
            F = griddedInterpolant({yI,xI},I(:,:,i),'linear','nearest');
            AI = F(Ys,Xs);
            AI(isnan(AI)) = 0; % should only be where W == 0
            
            % error
            err = AI - J{i};
            % gradient
            [AI_x,AI_y] = gradient2d(AI,dxJ(i,1),dxJ(i,2));
            Derr = zeros(size(AI,1),size(AI,2),6);
            count = 0;
            for r_ = 1 : 2
                for c = 1 : 3
                    dA = double((1:3==r_))' * double((1:3==c));
                    AdAAi = A*dA;
                    Xs = AdAAi(1,1)*XJ + AdAAi(1,2)*YJ + AdAAi(1,3);
                    Ys = AdAAi(2,1)*XJ + AdAAi(2,2)*YJ + AdAAi(2,3);
                    count = count + 1;
                    Derr(:,:,count) = AI_x.*Xs + AI_y.*Ys;
                end
            end
            tmp = reshape(bsxfun(@times,Derr,sqrt(W{i})),[],6);
            DerrDerr = tmp'*tmp;
            
            
            % update
            step = DerrDerr \ (tmp'*(err(:).*sqrt(W{i}(:))));
            step = reshape(step,3,2)';
            if any(isnan(step(:))) || any(isinf(step(:)))
                keyboard
                step = zeros(size(step));

            end
            % AJ = repmat(eye(3),[1,1,length(files)]);

            Ai(1:2,1:3) = Ai(1:2,1:3) - e * step;

            % rigid
            [U,S,V] = svd(Ai(1:2,1:2));
            Ai(1:2,1:2) = U*V';
            AJ(:,:,i) = inv(Ai);
        
            if i == round(size(I,3)/2)
                
                danfigure(25);
                clf
                subplot(2,2,1)
                imagesc(I(:,:,i));
                title('I')
                colorbar;
                
                subplot(2,2,2)
                imagesc(J{i})
                title('J');
                colorbar;
                
                subplot(2,2,3)
                imagesc(I0(:,:,i))
                title('I0')
                colorbar
                
                subplot(2,2,4)
                imagesc(err.*W{i})
                title('err')
                colorbar
                
                drawnow
            end
                
        
            end
        
        
        
        end
        
        if it == 1
            clim = quantile(I(:),qlim);
        end
        if ~(mod(it-1,5))
            danfigure(99);
            sliceView(xI,yI,zI,I,nplot,clim)
            danfigure(100);
            sliceView(xI,yI,zI,I0,nplot,clim)
            title(['down ' num2str(down) ', iter ' num2str(it)])

            % human dataset, produce output
            if contains(target_dir,'human')
                if it == 1 && downcount == 1
                    frame0 = [];
                    frame1 = [];
                end
                frame0 = [frame0, getframe(99)];
                frame1 = [frame1,getframe(100)];
                frame2Gif(frame0,[output_dir 'slice_average.gif'])
                frame2Gif(frame1,[output_dir 'slice_align.gif'])
            end

            drawnow;
        end
        % make it zero mean(?)
    end

    
    
end % end of down loop

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end
save([output_dir 'initializer_A.mat'],'AJ');

if contains(target_dir,'human')
    
    
    % write out volume    
    %write_vtk_image(xI,yI,zI,uint8(I0*255),[output_dir 'slice_recon.vtk'],'reconstructed human slices')
    addpath /home/dtward/Documents/Functions/avwQuiet
%     1    Binary             (  1 bit  per voxel)
%     2    Unsigned character (  8 bits per voxel)
%     4    Signed short       ( 16 bits per voxel)
%     8    Signed integer     ( 32 bits per voxel)
%    16    Floating point     ( 32 bits per voxel)
%    32    Complex, 2 floats  ( 64 bits per voxel), not supported
%    64    Double precision   ( 64 bits per voxel)
%   128    Red-Green-Blue     (128 bits per voxel), not supported
    
    avw = avw_hdr_make;
    avw.img = single(I0);
    avw.hdr.dime.dim(2:4) = size(I0);
    avw.hdr.dime.pixdim(2:4) = [dxJ(1,:),zJ0(2)-zJ0(1)]; 
    % note dx0 is arbitrary in this code, 
    % for human dz0 should be multiplied by 2 (for nissl/fluoro spacing),
    % or look at distance between nissl slices
    avw.hdr.dime.bitpix = 32;
    avw.hdr.dime.dtype = 16;
    avw.hdr.dime.bitpix = 8;
    avw.hdr.dime.dtype = 2;
    avw.img = uint8(I0*255);
    avw_img_write(avw,[output_dir 'slice_recon.img'])
    avw.img = uint8(I*255);
    avw_img_write(avw,[output_dir 'slice_average.img'])
    avw.img = uint8(W0*255);
    avw_img_write(avw,[output_dir 'slice_weights.img'])
end
% keyboard
return
keyboard


%%

for i = 1 : n
    danfigure(44);
    imagesc(I0(:,:,i))
    axis image
    title(['slice ' num2str(i) ' of ' num2str(n)])
    drawnow
%     pause
end