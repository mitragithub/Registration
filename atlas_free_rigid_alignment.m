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
for downcount = 1 : length(downs)
    down = downs(downcount);

    
    % loop through the files
    nxJ = zeros(length(files),2);
    for i = 1 : n
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
        J_(J_==1) = val(i);
        
        % for thick images, just ignore them
        if dz0(i) > skip_thick && skip_thick > 0
            if i > 1
                J_ = ones(size(J_))*val(i-1);
            else
                J_ = zeros(size(J_));
            end
        end
        
        
        % downsample
        [xJ{i},yJ{i},J{i}] = downsample2D(xJ{i},yJ{i},J_,down*[1 1]);
        % sometimes images will just be one pixel
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
        
        
    end
    nxI = [max(nxJ,[],1),length(files)];
    % last component doesn't matter, display as a square
    dxI = [dxJ(1,:),1/n*nxI(1)*dxJ0(1,1)*down]; % common to all slices
    I0 = zeros(nxI(2),nxI(1),nxI(3));
    

    
    xI = (0 : nxI(1)-1)*dxI(1);xI = xI - mean(xI);
    yI = (0 : nxI(2)-1)*dxI(2);yI = yI - mean(yI);
    zI = (0 : nxI(3)-1)*dxI(3);zI = zI - mean(zI);
    [XI_,YI_] = meshgrid(xI,yI);
    

    
    qlim = [0.001,0.999];
    clim = quantile(I0(:),qlim);
    
    %%
    % weights for averaging
    tmp = exp(-[-r:r].^2/2/(r/3)^2);
    tmp = tmp / sum(tmp);
%     CC = ([-r:r]==0) - tmp;
%     CC = CC*0.2;
    CC = tmp;
    CC(([-r:r]==0)) = 0;
    CC = CC/sum(CC);
    
    nmax = length(CC);
    
    
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
        end
        if it == 1
            qlim = [0.001,0.999];
            clim = quantile(I0(:),qlim);
        end

    
    
        % first pad the image
        npad = r;
        I0p = padarray(I0,[0,0,npad],'both','symmetric');
        I = convn(I0p,reshape(CC,1,1,[]),'same');
        I = I(:,:,npad+1:end-npad);
        % compute the image
        for i = 1 : n
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
            tmp = reshape(Derr,[],6);
            DerrDerr = tmp'*tmp;
            
            
            % update
            step = DerrDerr \ (tmp'*err(:));
            step = reshape(step,3,2)';
            if any(isnan(step(:))) || any(isinf(step(:)))
                step = zeros(size(step));
%                 keyboard
            end
            % AJ = repmat(eye(3),[1,1,length(files)]);

            Ai(1:2,1:3) = Ai(1:2,1:3) - e * step;

            % rigid
            [U,S,V] = svd(Ai(1:2,1:2));
            Ai(1:2,1:2) = U*V';
            AJ(:,:,i) = inv(Ai);
        
            if i == 1
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
                imagesc(err)
                title('err')
                colorbar
                
                drawnow
            end
                
        
        end
        
        danfigure(99);
        sliceView(xI,yI,zI,I)
        danfigure(100);
        sliceView(xI,yI,zI,I0)
        title(['down ' num2str(down) ', iter ' num2str(it)])
        
        
        
        drawnow;
        
        
        
        % make it zero mean(?)
    end

    
    
end % end of down loop

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end
save([output_dir 'initializer_A.mat'],'AJ');

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