function affine_for_initial_alignment(atlas_file, input_dir, pattern, detailed_output_dir, downs, niter)
% load slices
% construct 3D volume from initializer
% do affine, and update initializer


addpath Functions/plotting
addpath Functions/downsample
addpath Functions/vtk

% get the template
[xI,yI,zI,I,title_,names] = read_vtk_image(atlas_file);
I = double(I);


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

% % suppose this was wrong for DK39
% dxJ0 = dxJ0/4;

% reconstruct the slices
vars = load([detailed_output_dir 'initializer_A.mat']);
AJ = vars.AJ;

% initialize an image
% and one for averaging
if nargin < 5
downs = [8,4];
end
mindown = min(downs); % I can downsample it this much and put it in the next function
maxdown = max(downs); % I will need to pad it this much

if nargin < 6
    niter = 50;
end

%%
% loop over files and reconstruct slices into a 3D volume with approximatly
% 200 micron resolution (i.e. 50um allen atlas downsampled by 4)
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
        % typically we expect more slices than 50 micron spacing
        nslices = round((zJ0(end)-zJ0(1))/50/mindown);
        J = zeros([round(max(nxJ0(:,end:-1:1)/3/mindown)), nslices]);
        zJ = (0 : size(J,3)-1)*50*mindown;
        if nslices > length(files)
            nslices = length(files);
            J = zeros([round(max(nxJ0(:,end:-1:1)/3/mindown)), nslices]);            
            dz = (zJ0(end) - zJ0(1))/nslices;
            zJ = (0 : size(J,3)-1)*dz;
        end
        
        
        O = zeros(size(J));
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



% we need to pad this image before downsampling
% preprocess I, I want to pad it for example
dxI = [xI(2)-xI(1), yI(2)-yI(1), zI(2)-zI(1)];
for i = 1 : maxdown
I = padarray(I,[1,1,1],'both'); % zero pad
xI = [xI(1)-dxI(1), xI, xI(end)+dxI(1)];
yI = [yI(1)-dxI(2), yI, yI(end)+dxI(2)];
zI = [zI(1)-dxI(3), zI, zI(end)+dxI(3)];
end
[xI,yI,zI,I] = downsample(xI,yI,zI,I,mindown*[1,1,1]);





% scale it for numerical stability, since its scale doesn't matter
Imean = mean(I(:));
Istd = std(I(:));

I = I - mean(I(:));
I = I/std(I(:));

%%
% initial guess of affine transform
A = eye(4);
% permute xy, and flip z
% this is only valid for Allen vtk atlas
A = [0,1,0,0;
    1,0,0,0;
    0,0,-1,0;
    0,0,0,1]*A;
if contains(input_dir,'DK39')
    
%     A = eye(4);
% permute xy, and flip z
% this is only valid for Allen vtk atlas
% A = [0,1,0,0;
%     1,0,0,0;
%     0,0,-1,0;
%     0,0,0,1]*A;
A = eye(4);
% swap xz? no
% A = [0,0,1,0;
%     0,1,0,0;
%     1,0,0,0
%     0,0,0,1]*A; 
% sway xy? no
% A = [0,1,0,0;
%     1,0,0,0;
%     0,0,1,0
%     0,0,0,1]*A; 
% swap yz
A = [1,0,0,0;
    0,0,1,0;
    0,-1,0,0
    0,0,0,1]*A; 
% better, got the sagittal right
% might be good

end

%%

% better to start small here due to boundary issues
A = diag([0.8,0.8,0.8,1])*A;

% test the initial alignment
F = griddedInterpolant({yI,xI,zI},I,'linear','nearest');
[XJ,YJ,ZJ] = meshgrid(xJ,yJ,zJ);
Ai = inv(A);
Xs = Ai(1,1)*XJ + Ai(1,2)*YJ + Ai(1,3)*ZJ + Ai(1,4);
Ys = Ai(2,1)*XJ + Ai(2,2)*YJ + Ai(2,3)*ZJ + Ai(2,4);
Zs = Ai(3,1)*XJ + Ai(3,2)*YJ + Ai(3,3)*ZJ + Ai(3,4);
AI = F(Ys,Xs,Zs);
danfigure(5)
sliceView(xJ,yJ,zJ,J)
danfigure(6)
sliceView(xJ,yJ,zJ,AI)
% keyboard
% if nargin < 7
%     eT_factor = 2e-7;
% end
% if nargin < 8
%     eL_factor = 1e-14;
% end

%%
keyboard
A = ThreeD_to_3D_affine_registration_GN(xI,yI,zI,I,xJ,yJ,zJ,Jnorm,A,downs/mindown,niter,0.25);
save([detailed_output_dir 'initializer_A.mat'],'AJ','A');