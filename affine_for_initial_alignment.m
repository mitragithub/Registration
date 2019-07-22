function affine_for_initial_alignment(atlas_file, input_dir, pattern, detailed_output_dir)
% load slices
% construct 3D volume from initializer
% do affine, and update initializer

addpath Functions/plotting
addpath Functions/downsample
addpath Functions/vtk

% get the template
[xI,yI,zI,I,title_,names] = read_vtk_image(atlas_file);
I = double(I);
dxI = [xI(2)-xI(1), yI(2)-yI(1), zI(2)-zI(1)];
% scale it for numerical stability, since its scale doesn't matter
Imean = mean(I(:));
Istd = std(I(:));

I = I - mean(I(:));
I = I/std(I(:));



% get the target data
geometry_file = dir([input_dir '*.csv']);
fid = fopen([input_dir geometry_file(1).name],'rt');
line = fgetl(fid); % ignore the first line
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


% reconstruct the slices
vars = load([detailed_output_dir 'initializer_A.mat']);
AJ = vars.AJ;

% initialize an image
% and one for averaging

downs = [8,4];
mindown = min(downs);
[xI,yI,zI,I] = downsample(xI,yI,zI,I,mindown*[1,1,1]);

%%
% loop over files
for i = 1 : length(files)
    J_ = imread([input_dir files{i}]);
    J_ = mean(double(J_)/255.0,3);
    modeJ_ = mode(J_(J_(:)~=1));
    J_(J_==1) = modeJ_;
    xJ_ = (0 : nxJ0(i,1)-1)*dxJ0(i,1); xJ_ = xJ_ - mean(xJ_);
    yJ_ = (0 : nxJ0(i,2)-1)*dxJ0(i,1); yJ_ = yJ_ - mean(yJ_);
    
    % we want to downsample by 3 here to get resolution ~14 close to atlas ~45 
    [xJ_,yJ_,J_] = downsample2D(xJ_,yJ_,J_,[1,1]*3*mindown);
    
    danfigure(1);
    imagesc(xJ_,yJ_,J_)
    axis image;
    colormap gray;
    
    if i == 1
        J = zeros([round(max(nxJ0/3/mindown)), round((zJ0(end)-zJ0(1))/50/mindown)]);
        O = zeros(size(J));
        zJ = (0 : size(J,3)-1)*50*mindown;
        zJ = zJ - mean(zJ);
        xJ = (0 : size(J,2)-1)*(xJ_(2)-xJ_(1)); xJ = xJ - mean(xJ);
        yJ = (0 : size(J,1)-1)*(yJ_(2)-yJ_(1)); yJ = yJ - mean(yJ);
%         [XJ,YJ,ZJ] = meshgrid(xJ,yJ,zJ);
        [XJ_,YJ_] = meshgrid(xJ,yJ);
    end
    
    % apply transform
    F = griddedInterpolant({yJ_,xJ_},J_,'linear','none');
    B = AJ(:,:,i);
    Xs = B(1,1)*XJ_ + B(1,2)*YJ_ + B(1,3);
    Ys = B(2,1)*XJ_ + B(2,2)*YJ_ + B(2,3);
    J_R = F(Ys,Xs);
    J_R(isnan(J_R)) = modeJ_;
    danfigure(2);
    imagesc(xJ,yJ,J_R)
    axis image
    colormap gray;
    
    
    
    % find the slice this belongs to
    d2 = (zJ0(i) - zJ).^2;
    ind = find( d2==min(d2) ,1,'first');
    
    J(:,:,ind) = J(:,:,ind) + J_R;
    O(:,:,ind) = O(:,:,ind) + 1;
%     if ind > 1
%     J(:,:,ind-1) = J(:,:,ind-1) + J_R*0.5;
%     O(:,:,ind-1) = O(:,:,ind-1) + 0.5;
%     end
% 
%     if ind < size(J,3)
%     J(:,:,ind+1) = J(:,:,ind+1) + J_R*0.5;
%     O(:,:,ind+1) = O(:,:,ind+1) + 0.5;
%     end

    if ~mod(i-1,10) || i == length(files)
    danfigure(3);
    sliceView(xJ,yJ,zJ,J./O,5,[0,1])
%     danfigure(4);
%     sliceView(xJ,yJ,zJ,O,5,[0,1])
    
    
    drawnow
    end
    
    
    
end
Jnorm = J./O;
Jnorm(isnan(Jnorm)) = 0.8;
danfigure(3)
sliceView(xJ,yJ,zJ,Jnorm);
%%
close all
niter = 10;
A = eye(4);
% permute xy, and flip z
% this is only valid for Allen vtk atlas
A = [0,1,0,0;
    1,0,0,0;
    0,0,-1,0;
    0,0,0,1]*A;
eTfactor = 5e-7;
eLfactor = 2e-14;
A = ThreeD_to_3D_affine_registration(xI,yI,zI,I,xJ,yJ,zJ,Jnorm,A,downs/mindown,niter,eLfactor,eTfactor);
save([detailed_output_dir 'initializer_A.mat'],'AJ','A');