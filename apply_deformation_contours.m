function apply_deformation_contours(template_names,target_dir,detailed_output_dir,outdir)

disp(['Starting to apply deformations'])

% TODO: if transform is rigid, apply it using matrix, and instead of
% interpolation do translation, and rotation via shear
% NOTE: we apply a sequence of transformations to get from 3D atlas to 2D
% data
% x -> AJ . phiJ . A . phi . x
% originally we wanted to factor it as follows
% [AJ . phiJ] is registered to input
% [A . phi] is atlas to registered
% Instead, for each slice J we factored A into two parts
% say this introduces another linear transform B
% [AJ . phiJ . B] is registered to input
% [B^{-1} . A . phi] is atlas to registered
% B depends on the slice (linear in z)
%
% apply_deformation_contours repaces apply_deformatoin
% it computes grid lines and atlas labels 
% and it saves them in geo json format

% keyboard

addpath Functions/plotting
addpath Functions/vtk

rng('default') 
rng(1);
colors = rand(256,3);


if isempty(template_names)
    save_qc = 0;
else
    save_qc = 1;
end


%%
% load template
if save_qc
    disp('Starting to load templates')
    for c = 1 : length(template_names)
        disp(['Loaded template ' num2str(c) ' of ' num2str(length(template_names))])
        [xI{c},yI{c},zI{c},I{c},title_,names] = read_vtk_image(template_names{c});
        I{c} = double(I{c});
    end
    
end

%%
% make output directory
if ~exist(outdir,'dir')
    mkdir(outdir)
end

%%
% first thing is now to get slice thicknesses and location
disp('Starting to load slice geometry info')
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
    count = count + 1;
    % process this line, splitting at commas
    csv_data(count,:) = strsplit(line,',');
    %     
end
fclose(fid);
files = csv_data(:,1);
nxJ = cellfun(@(x)str2num(x), csv_data(:,2:3));
dxJ = cellfun(@(x)str2num(x), csv_data(:,5:6));

x0J = cellfun(@(x)str2num(x), csv_data(:,8:9));
z0J = cellfun(@(x)str2num(x), csv_data(:,10));

zJ = z0J;
dzJ = cellfun(@(x) str2num(x), csv_data(:,7));

for f = 1 : length(files)
    xJ{f} = x0J(f,1) + (0:nxJ(f,1)-1)*dxJ(f,1);
    yJ{f} = x0J(f,2) + (0:nxJ(f,2)-1)*dxJ(f,2);    
end
disp('finished loading slice geometry info')


%%
% get the transforms
% note that when working with 2 slice types, we'll have to load the
% combined transforms
disp('Starting to load saved data')
if exist([detailed_output_dir 'combined_A.mat'],'file')
Avars = load([detailed_output_dir 'combined_A.mat']);
vvars = load([detailed_output_dir 'combined_v.mat']);    
else
Avars = load([detailed_output_dir 'down1_A.mat']);
vvars = load([detailed_output_dir 'down1_v.mat']);
end
disp('Finished loading saved data')
A = Avars.A;
AJ = Avars.AJ;
dxVJ = zeros(length(files),2);
for f = 1 : length(files)
    dxVJ(f,:) = [vvars.xJ{f}(2)-vvars.xJ{f}(1),vvars.yJ{f}(2)-vvars.yJ{f}(1)];
end

vtx = vvars.vtx;
vty = vvars.vty;
vtz = vvars.vtz;

xV = vvars.xI;
yV = vvars.yI;
zV = vvars.zI;
[XV,YV,ZV] = meshgrid(xV,yV,zV);

nt = size(vtx,4);
dt = 1.0/nt;

vJtx = vvars.vJtx;
vJty = vvars.vJty;


%%
% now I have to flow the 3D transform
%% first the deformation
disp('Starting to integrate velocity field')
phiinvx = XV;
phiinvy = YV;
phiinvz = ZV;
for t = 1 : nt
    disp(['Integrating v field step '  num2str(t) ' of ' num2str(nt)])
    % update phi
    Xs = XV - vtx(:,:,:,t)*dt;
    Ys = YV - vty(:,:,:,t)*dt;
    Zs = ZV - vtz(:,:,:,t)*dt;
    % subtract and add identity
    F = griddedInterpolant({yV,xV,zV},phiinvx-XV,'linear','nearest');
    phiinvx = F(Ys,Xs,Zs) + Xs;
    F = griddedInterpolant({yV,xV,zV},phiinvy-YV,'linear','nearest');
    phiinvy = F(Ys,Xs,Zs) + Ys;
    F = griddedInterpolant({yV,xV,zV},phiinvz-ZV,'linear','nearest');
    phiinvz = F(Ys,Xs,Zs) + Zs;
end

phix = XV;
phiy = YV;
phiz = ZV;
disp('Starting to integrate velocity field backward')
for t = nt : -1 : 1
    disp(['Integrating v field step '  num2str(t) ' of ' num2str(nt)])
    % update phi
    Xs = XV + vtx(:,:,:,t)*dt;
    Ys = YV + vty(:,:,:,t)*dt;
    Zs = ZV + vtz(:,:,:,t)*dt;
    % subtract and add identity
    F = griddedInterpolant({yV,xV,zV},phix-XV,'linear','nearest');
    phix = F(Ys,Xs,Zs) + Xs;
    F = griddedInterpolant({yV,xV,zV},phiy-YV,'linear','nearest');
    phiy = F(Ys,Xs,Zs) + Ys;
    F = griddedInterpolant({yV,xV,zV},phiz-ZV,'linear','nearest');
    phiz = F(Ys,Xs,Zs) + Zs;
end
Aphix = A(1,1)*phix + A(1,2)*phiy + A(1,3)*phiz + A(1,4);
Aphiy = A(2,1)*phix + A(2,2)*phiy + A(2,3)*phiz + A(2,4);
Aphiz = A(3,1)*phix + A(3,2)*phiy + A(3,3)*phiz + A(3,4);
% get det jac
dxV = [xV(2)-xV(1), yV(2)-yV(1), zV(2)-zV(1)];
[phix_x,phix_y,phix_z] = gradient(phix,dxV(1),dxV(2),dxV(3));
[phiy_x,phiy_y,phiy_z] = gradient(phiy,dxV(1),dxV(2),dxV(3));
[phiz_x,phiz_y,phiz_z] = gradient(phiz,dxV(1),dxV(2),dxV(3));
detjac = (phix_x.*(phiy_y.*phiz_z - phiy_z.*phiz_y) ...
    - phix_y.*(phiy_x.*phiz_z - phiy_z.*phiz_x) ...
    + phix_z.*(phiy_x.*phiz_y - phiy_y.*phiz_x))*det(A);
% write these out
disp(['Writing out saved transformations and jacobians'])
if ~exist(outdir,'dir');mkdir(outdir);end;
% don't use this data for atlas to registered, the code below incorporates
% corrections that center the data for display in the MBA portal
% write_vtk_image(xV,yV,zV,single(cat(4,Aphix-XV,Aphiy-YV,Aphiz-ZV)),[outdir 'atlas_to_registered_displacement.vtk'],'atlas_to_registered')
write_vtk_image(xV,yV,zV,single(detjac),[outdir 'atlas_to_registered_detjac.vtk'],'atlas_to_registered_detjac')
% note this additional transform is a shear so does not impact detjac

% we also want the velocity field
write_vtk_image(xV,yV,zV,single(permute(cat(5,vtx,vty,vtz),[1,2,3,5,4])),[outdir 'atlas_to_registered_velocity.vtk'],'atlas_to_registered_velocity')
disp(['Finished writing out saved transformations and jacobians'])

%%
% this only makes sense when labels and deformation are at the same size
% and resolution
% detjacbrain = detjac((padarray(I{1},[1,1,1])>0)).^(1/3);
% meandetjacbrain = mean(detjacbrain);
% disp(['Mean change of brain scale is ' num2str(meandetjacbrain)])
% mean_absolute_deviation = mean(abs(detjacbrain - meandetjacbrain));
% disp(['MAD is ' num2str(mean_absolute_deviation)])
% return


%%
% update atlas-to-registered based on in-plane shifts that were calculated
% to center data for display
% this section corrects for angle, and xz and yz shear

% unit vector up
v = [1;0;0;0];
% what happens to unit vector
Av = A*v;
theta = atan2(-Av(1),Av(2));
% what rotation corresponded to this?
R = [cos(theta),-sin(theta),0,0;
    sin(theta),cos(theta),0,0;
    0,0,1,0;
    0,0,0,1];
% and we also have a shift that depended on z linearly
%   
u0 = [0;0;zJ(f);1]; % what slice am I on?
v0 = A\u0; % where does this appear in the atlas
zatlas = v0(3); % what z location does it appear?
% look at Pz
Pz = [0,0,0,0;
    0,0,0,0;
    0,0,1,0;
    0,0,0,1];
Pxy = [1,0,0,0;
    0,1,0,0;
    0,0,0,0;
    0,0,0,1];
% the goal here is to add a z depending xy shift
SHIFT = (Pxy * A * Pz * inv(A) * Pz);
SHIFT(1,1) = 1;
SHIFT(2,2) = 1;
SHIFT(3,3) = 1;

% update A and AJ
A_ = inv(R)*inv(SHIFT)*A;
AJ_ = zeros(size(AJ));
for i = 1 : size(AJ,3)
    SHIFTz = eye(3);
    SHIFTz(1,3) = SHIFT(1,4) + SHIFT(1,3)*zJ(i);
    SHIFTz(2,3) = SHIFT(2,4) + SHIFT(2,3)*zJ(i);
    AJ_(:,:,i) = AJ(:,:,i) * SHIFTz * R([1,2,4],[1,2,4]) ;
end

A_phix = A_(1,1)*phix + A_(1,2)*phiy + A_(1,3)*phiz + A_(1,4);
A_phiy = A_(2,1)*phix + A_(2,2)*phiy + A_(2,3)*phiz + A_(2,4);
A_phiz = A_(3,1)*phix + A_(3,2)*phiy + A_(3,3)*phiz + A_(3,4);
write_vtk_image(xV,yV,zV,single(cat(4,A_phix-XV,A_phiy-YV,A_phiz-ZV)),[outdir 'atlas_to_registered_displacement.vtk'],'atlas_to_registered')


%%
% now I have to loop through slices
% for f = 681%1 : 1 : length(files)
for f = 1 : 1 : length(files)
    disp(['Starting to apply transforms for slice ' num2str(f) ' of ' num2str(length(files))]);
    [dir_,fname_,ext_] = fileparts(files{f});
    % when writing vtk, I'll calculate dx from the "z" variable
    % so we need 2 elements
    zJ_write = zJ(f) + [0,1]*dzJ(f);
    xVJ = vvars.xJ{f};
    yVJ = vvars.yJ{f};
    [XVJ,YVJ] = meshgrid(xVJ,yVJ);
    
    %%
    % now I have to flow the 2D transform
    ntJ = size(vJtx{f},3);
    dtJ = 1/ntJ;
    
    phiJinvx = XVJ;
    phiJinvy = YVJ;
    for t = 1 : ntJ
        % update phi
        Xs = XVJ - vJtx{f}(:,:,t)*dtJ;
        Ys = YVJ - vJty{f}(:,:,t)*dtJ;
        % subtract and add identity
        F = griddedInterpolant({yVJ,xVJ},phiJinvx-XVJ,'linear','nearest');
        phiJinvx = F(Ys,Xs) + Xs;
        F = griddedInterpolant({yVJ,xVJ},phiJinvy-YVJ,'linear','nearest');
        phiJinvy = F(Ys,Xs) + Ys;
    end
    
    %%
    % now I have to get coordinates for this image
    % we will apply the transformations to these points
    J = imread([target_dir files{f}]);
    [XJ,YJ] = meshgrid(xJ{f},yJ{f});
%     danfigure(1);
%     imagesc(xJ{f},yJ{f},J)
%     axis image


    %%
    % first apply the 2D affine
    BJ = inv(AJ(:,:,f));
    AJiX = BJ(1,1)*XJ + BJ(1,2)*YJ + BJ(1,3);
    AJiY = BJ(2,1)*XJ + BJ(2,2)*YJ + BJ(2,3);
    AJiZ = zJ(f); % identity
    
    %%
    % now apply 2D diffeo
    F = griddedInterpolant({yVJ,xVJ}, phiJinvx-XVJ, 'linear','nearest');
    phiJiAJiX = F(AJiY,AJiX) + AJiX;
    F = griddedInterpolant({yVJ,xVJ}, phiJinvy-YVJ, 'linear','nearest');
    phiJiAJiY = F(AJiY,AJiX) + AJiY;
    phiJiAJiZ = AJiZ; % identity
    
    %%
    % now apply the affine
    B = inv(A);
    AiPhiJiAJiX = B(1,1)*phiJiAJiX + B(1,2)*phiJiAJiY + B(1,3)*phiJiAJiZ + B(1,4);
    AiPhiJiAJiY = B(2,1)*phiJiAJiX + B(2,2)*phiJiAJiY + B(2,3)*phiJiAJiZ + B(2,4);
    AiPhiJiAJiZ = B(3,1)*phiJiAJiX + B(3,2)*phiJiAJiY + B(3,3)*phiJiAJiZ + B(3,4);
    
    %%
    % now apply the diffeo
    F = griddedInterpolant({yV,xV,zV}, phiinvx-XV, 'linear','nearest');
    phiiAiPhiJiAJiX = F(AiPhiJiAJiY,AiPhiJiAJiX,AiPhiJiAJiZ) + AiPhiJiAJiX;
    F = griddedInterpolant({yV,xV,zV}, phiinvy-YV, 'linear','nearest');
    phiiAiPhiJiAJiY = F(AiPhiJiAJiY,AiPhiJiAJiX,AiPhiJiAJiZ) + AiPhiJiAJiY;
    F = griddedInterpolant({yV,xV,zV}, phiinvz-ZV, 'linear','nearest');
    phiiAiPhiJiAJiZ = F(AiPhiJiAJiY,AiPhiJiAJiX,AiPhiJiAJiZ) + AiPhiJiAJiZ;
    
    
    %%
    % now I have to apply the transform to my atlas
    % we don't need to make this figure, just the "straight" version
    
    
    %%
    % we also need the "reconstruction" coordinates (backwards for slices)
    % make it the same size as the full res image
    phi1tJinvx = XVJ;
    phi1tJinvy = YVJ;
    for t = ntJ : -1 : 1 % I think I had this order wrong last time!
        % update phi
        Xs = XVJ + vJtx{f}(:,:,t)*dtJ;
        Ys = YVJ + vJty{f}(:,:,t)*dtJ;
        % subtract and add identity
        F = griddedInterpolant({yVJ,xVJ},phi1tJinvx-XVJ,'linear','nearest');
        phi1tJinvx = F(Ys,Xs) + Xs;
        F = griddedInterpolant({yVJ,xVJ},phi1tJinvy-YVJ,'linear','nearest');
        phi1tJinvy = F(Ys,Xs) + Ys;
    end
    % the reconstruction is gonna be 
    % AphiJ
    AJphiJX = AJ(1,1)*phi1tJinvx + AJ(1,2)*phi1tJinvy + AJ(1,3);
    AJphiJY = AJ(2,1)*phi1tJinvx + AJ(2,2)*phi1tJinvy + AJ(2,3);
    % upsample onto size of full res image
    F = griddedInterpolant({yVJ,xVJ},AJphiJX-XVJ,'linear','nearest');
    AJphiJX = F(YJ,XJ) + XJ;
    F = griddedInterpolant({yVJ,xVJ},AJphiJY-YVJ,'linear','nearest');
    AJphiJY = F(YJ,XJ) + YJ;
    
    % I don't need to save these because I'm only saving the straightened
    % version

    
    %%
    % factor the matrix A to make sure it is z is pointing up and there is
    % no shear
    u0 = [0;0;zJ(f);1];
    v0 = A\u0;
    zatlas = v0(3);
    tmp = A*[0;0;zatlas;1];
    xyoff = tmp(1:2);
    
    % unit vector up
    v = [1;0;0;0];
    % what happens to unit vector
    Av = A*v;
    
    % so we want B or Axy to do two things
    % correct this translation (in xy)
    % and create this rotation (y component)
    theta = atan2(-Av(1),Av(2));
    % what rotation corresponded to this?
    R = [cos(theta),-sin(theta),0,0;
        sin(theta),cos(theta),0,0;
        0,0,1,0;
        0,0,0,1];
    
%     inv(R)*A*v has x component 0
    
    % first apply A
    % then shift back to center
    % then rotate to verticle
    % then inverse rotate
    % then inverse shift
    Shift = [1,0,0,xyoff(1);
        0,1,0,xyoff(2);
        0,0,1,0;
        0,0,0,1];
    Axyz = inv(R)*inv(Shift)*A;
    Axy = Shift*R;
    
    % get a 2D version    
    Axy_ = Axy([1,2,4],[1,2,4]);
    

    % so the better reconstruction
    % x \to Axy x \to phiJ(Axy(x)) \to AJ(phiJ(Axy(x)));
    % we need this forward transform to pull the image J backwards
    AxyX = Axy_(1,1)*XJ + Axy_(1,2)*YJ + Axy_(1,3);
    AxyY = Axy_(2,1)*XJ + Axy_(2,2)*YJ + Axy_(2,3);
    % now evaluate phiJ here
    F = griddedInterpolant({yVJ,xVJ},phi1tJinvx-XVJ,'linear','nearest');
    phiJAxyX = F(AxyY,AxyX) + AxyX;
    F = griddedInterpolant({yVJ,xVJ},phi1tJinvy-YVJ,'linear','nearest');
    phiJAxyY = F(AxyY,AxyX) + AxyY;
    % now apply AJ
    AJphiJAxyX = AJ(1,1,f)*phiJAxyX + AJ(1,2,f)*phiJAxyY + AJ(1,3,f);
    AJphiJAxyY = AJ(2,1,f)*phiJAxyX + AJ(2,2,f)*phiJAxyY + AJ(2,3,f);
    % let's deform J with this as a test
    Jrecon = zeros(size(J));
    for c = 1 : size(J,3)
        F = griddedInterpolant({yJ{f},xJ{f}},double(J(:,:,c))/255.0,'linear','nearest');
        Jrecon(:,:,c) = F(AJphiJAxyY,AJphiJAxyX);
    end


    % the transform to calculate Jrecon is a map from registered to input
    out = single(cat(4,AJphiJAxyX-XJ,AJphiJAxyY-YJ, zeros(size(AJphiJAxyX))));
    write_vtk_image(xJ{f},yJ{f},zJ_write,out,[outdir 'registered_to_input_displacement_' fname_ '.vtk'],'registered_to_input_displacement');

    
%     %%
    % now we need a transform from input to registered, its the inverse of
    % the above
    % y = AJ(phiJ(Axy(x)))
    % x = Axyi(phiJi(AJi(y)))
    % i have  phiJiAJiX 
    % just need to apply Axyi
    AxyiphiJiAjiX = Axy_(1,1)*phiJiAJiX + Axy_(1,2)*phiJiAJiY + Axy_(1,3);
    AxyiphiJiAjiY = Axy_(2,1)*phiJiAJiX + Axy_(2,2)*phiJiAJiY + Axy_(2,3);
    out = single(cat(4,AxyiphiJiAjiX-XJ,AxyiphiJiAjiY-YJ, zeros(size(AxyiphiJiAjiY))));
    write_vtk_image(xJ{f},yJ{f},zJ_write,out,[outdir 'input_to_registered_displacement_' fname_ '.vtk'],'input_to_registered_displacement');

    
    
%     %%
    % and now we need the deformation
    % do it exactly the same as above, but use identity
    BJ = eye(3);
    AJiX = BJ(1,1)*XJ + BJ(1,2)*YJ + BJ(1,3);
    AJiY = BJ(2,1)*XJ + BJ(2,2)*YJ + BJ(2,3);
    AJiZ = zJ(f); % identity
    
    % now apply 2D diffeo
    F = griddedInterpolant({yVJ,xVJ}, (phiJinvx-XVJ)*0, 'linear','nearest');
    phiJiAJiX = F(AJiY,AJiX) + AJiX;
    F = griddedInterpolant({yVJ,xVJ}, (phiJinvy-YVJ)*0, 'linear','nearest');
    phiJiAJiY = F(AJiY,AJiX) + AJiY;
    phiJiAJiZ = AJiZ; % identity
    
    % now apply the affine
    % it should be called Ai, not Bi
    B = inv(Axyz);
    AxyziPhiJiAJiX = B(1,1)*phiJiAJiX + B(1,2)*phiJiAJiY + B(1,3)*phiJiAJiZ + B(1,4);
    AxyziPhiJiAJiY = B(2,1)*phiJiAJiX + B(2,2)*phiJiAJiY + B(2,3)*phiJiAJiZ + B(2,4);
    AxyziPhiJiAJiZ = B(3,1)*phiJiAJiX + B(3,2)*phiJiAJiY + B(3,3)*phiJiAJiZ + B(3,4);
    
    % now apply the diffeo
    F = griddedInterpolant({yV,xV,zV}, phiinvx-XV, 'linear','nearest');
    phiiAxyziPhiJiAJiX = F(AxyziPhiJiAJiY,AxyziPhiJiAJiX,AxyziPhiJiAJiZ) + AxyziPhiJiAJiX;
    F = griddedInterpolant({yV,xV,zV}, phiinvy-YV, 'linear','nearest');
    phiiAxyziPhiJiAJiY = F(AxyziPhiJiAJiY,AxyziPhiJiAJiX,AxyziPhiJiAJiZ) + AxyziPhiJiAJiY;
    F = griddedInterpolant({yV,xV,zV}, phiinvz-ZV, 'linear','nearest');
    phiiAxyziPhiJiAJiZ = F(AxyziPhiJiAJiY,AxyziPhiJiAJiX,AxyziPhiJiAJiZ) + AxyziPhiJiAJiZ;    
    

    
    
%     %%
    % now we want a deformation that includes straightening
    qlim = [0.01,0.99];

    if save_qc
        for t = 1 : length(template_names)
            % if the template is an allen annotation, we do it with
            % outlines
            
            % if the template is waxholm annotation, we do it with outlines
            danfigure(5);
            clf;
            [thedir_,thename_,theext_] = fileparts(template_names{t});
            if ~isempty(strfind(template_names{t},'annotation')) || ~isempty(strfind(thename_,'rat_atlas'))
                
                F = griddedInterpolant({yI{t},xI{t},zI{t}},I{t},'nearest','nearest');
                Seg = F(phiiAxyziPhiJiAJiY, phiiAxyziPhiJiAJiX, phiiAxyziPhiJiAJiZ);
%                 [Seg_x,Seg_y] = gradient(Seg);
%                 Seg_contour = Seg_x~=0 | Seg_y~=0;
%                 % dilate it
%                 Seg_contour_d = (convn(Seg_contour,ones(5),'same')~=0) - Seg_contour;
% 
%                 
%                 %             SegInd = mod(Seg,256)+1;
%                 %             SegRGB = cat(3,reshape(colors(SegInd,1),size(Seg)),...
%                 %                 reshape(colors(SegInd,2),size(Seg)),...
%                 %                 reshape(colors(SegInd,3),size(Seg)));
%                 %             imagesc(xJ{f},yJ{f},SegRGB*segalpha + Jrecon*(1-segalpha))
%                 % plot the contour instead
%                 Jshow = Jrecon;
%                 tmp = mean(Jshow,3);
%                 clim = quantile(tmp(:),qlim);
%                 Jshow = (Jshow - clim(1))/diff(clim);
%                 for c = 1 : 3
%                     tmp = Jshow(:,:,c);
%                     tmp(Seg_contour~=0) = 1;
%                     tmp(Seg_contour_d~=0) = 0;
%                     Jshow(:,:,c) = tmp;
%                 end
%                 imagesc(xJ{f},yJ{f},Jshow);axis image

                
                
                
                % get the unique labels
 
                % this isn't any good
                % we need to blur the labels to their native resolution
                res0 = xI{t}(2)-xI{t}(1);
                res1 = xJ{f}(2)-xJ{f}(1);
                pixblur = res0/res1;
                rblur = ceil(pixblur/2);
                rblur = ceil(pixblur);
                % stack up the images with shifts
                Seg_ = [];
                for r = -rblur : rblur
                    for c = -rblur : rblur
                        if r^2 + c^2 > rblur^2
                            continue
                        end
                        Seg_(:,:,end+1) = circshift(Seg,[r,c]);
                    end
                end
                % take the mode
                out = mode(Seg_,3);
                Seg = out;
                      
                [Seg_x,Seg_y] = gradient(Seg);
                Seg_contour = Seg_x~=0 | Seg_y~=0;
                % dilate it
                Seg_contour_d = (convn(Seg_contour,ones(5),'same')~=0) - Seg_contour;

                
                Jshow = Jrecon;
                tmp = mean(Jshow,3);
                clim = quantile(tmp(:),qlim);
                Jshow = (Jshow - clim(1))/diff(clim);
                for c = 1 : 3
                    tmp = Jshow(:,:,c);
                    tmp(Seg_contour~=0) = 1;
                    tmp(Seg_contour_d~=0) = 0;
                    Jshow(:,:,c) = tmp;
                end
%                 imagesc(xJ{f},yJ{f},Jshow);axis image
                imagesc(xJ{f},yJ{f},Jrecon);axis image
                
                
%                 structure_contours = {};
%                 labels = unique(Seg(:));
%                 for l = 1 : length(labels)
%                     out = contour(xJ{f},yJ{f},Seg==labels(l));
%                 end
                
                % let's use the contour approach and do it a few ways
                % method 1, previous version
                % define the interpolator
                
                % resample to 0.3 um
%                 keyboard
%                 scale1 = 2.0;
%                 x1 = xJ{f}(1):scale1:xJ{f}(end);
%                 y1 = yJ{f}(1):scale1:yJ{f}(end);
%                 [X1,Y1] = meshgrid(x1,y1);
%                 F = griddedInterpolant({yJ{f},xJ{f}},phiiAxyziPhiJiAJiX,'nearest','nearest');
%                 phiiX1 = F(Y1,X1);
%                 F = griddedInterpolant({yJ{f},xJ{f}},phiiAxyziPhiJiAJiY,'nearest','nearest');
%                 phiiY1 = F(Y1,X1);
%                 F = griddedInterpolant({yJ{f},xJ{f}},phiiAxyziPhiJiAJiZ,'nearest','nearest');
%                 phiiZ1 = F(Y1,X1);
%                 F = griddedInterpolant({yI{t},xI{t},zI{t}},I{t},'nearest','nearest');
%                 Seg = F(phiiY1, phiiX1, phiiZ1);
%                 imagesc(x1,y1,mod(Seg,64))
%                 structure_contours = {};                
%                 labels = unique(Seg(:));
%                 for l = 1 : length(labels)
%                     structure_contours{l} = contourc(x1,y1,double(Seg==labels(l)),[0.5,0.5]); % draw single level contours
%                 end
%                 close all;
%                 danfigure(999)
%                 clf;
%                 hold on;
%                 % plot them
%                 for l = 1 : length(labels)
%                     sc = structure_contours{l};
%                     while 1
%                         if isempty(sc)
%                             break
%                         end
%                         n = sc(2,1);
%                         xdata = sc(1,2:n+1);
%                         ydata = sc(2,2:n+1);
%                         sc = sc(:,n+2:end);
%                         plot(xdata,ydata)
%                     end
%                 end
%                 axis image
%                 title('Contours after high res resampling')
%                 set(gca,'ydir','reverse')
%                 saveas(999,'method1_out.png')
%                 set(gca,'xlim',[0 1000],'ylim',[0 1000])
%                 saveas(999,'method1_in.png')
%                 
%                 
%                 
%                 % method 2 is to draw contours after low resolution
%                 % resampling
%                 scale1 = 50; % I"m using 50 micron segs
%                 x1 = xJ{f}(1):scale1:xJ{f}(end);
%                 y1 = yJ{f}(1):scale1:yJ{f}(end);
%                 [X1,Y1] = meshgrid(x1,y1);
%                 F = griddedInterpolant({yJ{f},xJ{f}},phiiAxyziPhiJiAJiX,'nearest','nearest');
%                 phiiX1 = F(Y1,X1);
%                 F = griddedInterpolant({yJ{f},xJ{f}},phiiAxyziPhiJiAJiY,'nearest','nearest');
%                 phiiY1 = F(Y1,X1);
%                 F = griddedInterpolant({yJ{f},xJ{f}},phiiAxyziPhiJiAJiZ,'nearest','nearest');
%                 phiiZ1 = F(Y1,X1);
%                 F = griddedInterpolant({yI{t},xI{t},zI{t}},I{t},'nearest','nearest');
%                 Seg = F(phiiY1, phiiX1, phiiZ1);
% %                 imagesc(x1,y1,mod(Seg,64))
%                 structure_contours = {};                
%                 labels = unique(Seg(:));
%                 for l = 1 : length(labels)
%                     structure_contours{l} = contourc(x1,y1,double(Seg==labels(l)),[0.5,0.5]); % draw single level contours
%                 end
%                 close all;
%                 danfigure(999)
%                 clf
%                 hold on;
%                 % plot them
%                 for l = 1 : length(labels)
%                     sc = structure_contours{l};
%                     while 1
%                         if isempty(sc)
%                             break
%                         end
%                         n = sc(2,1);
%                         xdata = sc(1,2:n+1);
%                         ydata = sc(2,2:n+1);
%                         sc = sc(:,n+2:end);
%                         plot(xdata,ydata)
%                     end
%                 end
%                 axis image
%                 title('Contours after Low res resampling')
%                 set(gca,'ydir','reverse')
%                 saveas(999,'method2_out.png')
%                 set(gca,'xlim',[0 1000],'ylim',[0 1000])
%                 saveas(999,'method2_in.png')
% 
%                 
%                 % method 3 is to use high res
%                 % but use special interpolation
%                 % no need to do higher than 10
%                 scale1 = 10.0;
%                 x1 = xJ{f}(1):scale1:xJ{f}(end);
%                 y1 = yJ{f}(1):scale1:yJ{f}(end);
%                 [X1,Y1] = meshgrid(x1,y1);
%                 F = griddedInterpolant({yJ{f},xJ{f}},phiiAxyziPhiJiAJiX,'nearest','nearest');
%                 phiiX1 = F(Y1,X1);
%                 F = griddedInterpolant({yJ{f},xJ{f}},phiiAxyziPhiJiAJiY,'nearest','nearest');
%                 phiiY1 = F(Y1,X1);
%                 F = griddedInterpolant({yJ{f},xJ{f}},phiiAxyziPhiJiAJiZ,'nearest','nearest');
%                 phiiZ1 = F(Y1,X1);
%                 
%                 
% %                 F = griddedInterpolant({yI{t},xI{t},zI{t}},I{t},'nearest','nearest');
% %                 Seg = F(phiiY1, phiiX1, phiiZ1);
%                 if ~exist('unique_labels','var')
%                     disp('finding unique labels')
%                     unique_labels = unique(I{t}(:)); % this should only be done once!
%                 end
%                 Seg = zeros(size(phiiY1));
%                 SegP = zeros(size(phiiY1));
%                 for u = 1 : length(unique_labels)
%                     disp(['label ' num2str(u) ' of ' num2str(length(unique_labels))])
%                     F = griddedInterpolant({yI{t},xI{t},zI{t}},double(I{t}==unique_labels(u)),'nearest','nearest');
%                     SegP_ = F(phiiY1,phiiX1,phiiZ1);
%                     ind = SegP_>SegP;
%                     Seg(ind) = unique_labels(u);
%                     SegP = max(SegP,SegP_);
%                     imagesc(x1,y1,mod(Seg,64))
%                     drawnow;
%                 end
%                 
%                 
%                 
%                 imagesc(x1,y1,mod(Seg,64))
%                 structure_contours = {};                
%                 labels = unique(Seg(:));
%                 for l = 1 : length(labels)
%                     structure_contours{l} = contourc(x1,y1,double(Seg==labels(l)),[0.5,0.5]); % draw single level contours
%                 end
%                 close all;
%                 danfigure(999)
%                 clf;
%                 hold on;
%                 % plot them
%                 for l = 1 : length(labels)
%                     sc = structure_contours{l};
%                     while 1
%                         if isempty(sc)
%                             break
%                         end
%                         n = sc(2,1);
%                         xdata = sc(1,2:n+1);
%                         ydata = sc(2,2:n+1);
%                         sc = sc(:,n+2:end);
%                         plot(xdata,ydata)
%                     end
%                 end
%                 axis image
%                 title('Contours after fancy resampling')
%                 set(gca,'ydir','reverse')
%                 saveas(999,'method3_out.png')
%                 set(gca,'xlim',[0 1000],'ylim',[0 1000])
%                 saveas(999,'method3_in.png')
                
                
                
                % I could also blur but I want to smooth the points
                
                % method 4 is to draw contours after normal sampling
%                 keyboard
                % and smooth them                
                % oversample just a bit
                scale1 = 5; % I"m using 50 micron segs
                x1 = xJ{f}(1):scale1:xJ{f}(end);
                y1 = yJ{f}(1):scale1:yJ{f}(end);
                [X1,Y1] = meshgrid(x1,y1);
                F = griddedInterpolant({yJ{f},xJ{f}},phiiAxyziPhiJiAJiX,'nearest','nearest');
                phiiX1 = F(Y1,X1);
                F = griddedInterpolant({yJ{f},xJ{f}},phiiAxyziPhiJiAJiY,'nearest','nearest');
                phiiY1 = F(Y1,X1);
                F = griddedInterpolant({yJ{f},xJ{f}},phiiAxyziPhiJiAJiZ,'nearest','nearest');
                phiiZ1 = F(Y1,X1);
                F = griddedInterpolant({yI{t},xI{t},zI{t}},I{t},'nearest','nearest');
                Seg = F(phiiY1, phiiX1, phiiZ1);
%                 imagesc(x1,y1,mod(Seg,64))
                structure_contours = {};                
                labels = unique(Seg(:));
                for l = 1 : length(labels)
                    structure_contours{l} = contourc(x1,y1,double(Seg==labels(l)),[0.5,0.5]); % draw single level contours
                end
%                 close all;
%                 danfigure(999)
%                 clf
                hold on;
                % plot them
                ids = [];
                contours = {};
                for l = 1 : length(labels)
                    ids(l) = labels(l);
                    sc = structure_contours{l};
                    count = 0;
                    while 1
                        if isempty(sc)
                            break
                        end
                        count = count + 1;
                        n = sc(2,1);
                        xdata = sc(1,2:n+1);
                        ydata = sc(2,2:n+1);
                        
%                         hold off; 
%                         plot(xdata,ydata);
%                         hold on;
                        
                        % note the following
                        % x data and ydata have the same first and last
                        % point
                        % I can smooth things by projecting onto a family
                        % of functions
                        % just consider the x component and the y component
                        s = 1 : (length(xdata)-1);
                        s = s - mean(s);
                        b = exp(-s.^2/2.0/20^2);
                        b = b / sum(b.^2);
%                         B = zeros(length(s),length(s));
%                         for bbb = 1 : length(s)
%                             B(:,bbb) = circshift(b,bbb);
%                         end
                        v = b;
                        B = toeplitz([v(1) fliplr(v(2:end))], v);
%                         smoothdown = 5;
%                         B = B(:,1:smoothdown:end);
                        smooth = 1e-3;
                        BB = transpose(B)*B + eye(size(B,2)) * smooth;
%                         coeffsx = BB \ (transpose(B)* transpose(xdata(1:end-1)));
%                         xdatanew = B*coeffsx;
%                         coeffsy = BB \ (transpose(B)* transpose(ydata(1:end-1)));                        
%                         ydatanew = B*coeffsy;
                        coeffs = BB \ (transpose(B) * [xdata(1:end-1);ydata(1:end-1)]');
                        xdatanew = B*coeffs(:,1);
                        ydatanew = B*coeffs(:,2);
                        
                        % theoptimization
                        % c = arg_min |c_i b(x - x_i) - f(x)|^2 
                        % plus some regularization
                        
                        % or
%                         bhat = fft(b);
%                         xdatanew_ = ifft( (fft(xdata(1:end-1)) .* bhat) ./ (abs(bhat).^2 + 0.1 ) )';
%                         ydatanew_ = ifft( (fft(ydata(1:end-1)) .* bhat) ./ (abs(bhat).^2 + 0.1 ) )';
%                         plot([xdatanew_; xdatanew_(1)],[ydatanew_; ydatanew_(1)])


                        sc = sc(:,n+2:end);
                        
                        
                        
                        
                        % if it is too small, do not add it
                        if length(xdatanew) < 20
                            count = count - 1;
                            continue
                        end
                        plot([xdatanew; xdatanew(1)],[ydatanew; ydatanew(1)],'k','linewidth',2)
                        plot([xdatanew; xdatanew(1)],[ydatanew; ydatanew(1)],'w')
%                         drawnow
                        
                        
                        contours{l}{count} = [[xdatanew; xdatanew(1)],[ydatanew; ydatanew(1)]]';
                        
                        
                        
                    end
                    
                end
%                 axis image
%                 title('Contours after Low res resampling and splines')
%                 set(gca,'ydir','reverse')
%                 saveas(999,'method4_out.png')
%                 set(gca,'xlim',[0 1000],'ylim',[0 1000])
%                 saveas(999,'method4_in.png')
                
                
                
                % write out contours in geojson
                [directory_, filename_, extension_] = fileparts(files{f});
                contourname = [outdir filename_  '_contour_straight.geojson'];
                
                % note it is possible that ids will be longer than contours
                % if the contours are empty
                write_geojson_ontology(contourname,ids(1:length(contours)),contours);
                
                
  
            else
                
                F = griddedInterpolant({yI{t},xI{t},zI{t}},I{t},'linear','nearest');
                Seg = F(phiiAxyziPhiJiAJiY, phiiAxyziPhiJiAJiX, phiiAxyziPhiJiAJiZ);
                
                
                % otherwise, if not an annotation, we do the following
%                 imagesc(xJ{f},yJ{f},bsxfun(@plus,Seg/max(Seg(:))*segalpha, Jrecon*(1-segalpha)))
                
                
%                 J1 = mean(Jrecon,3);
%                 J1 = J1 - min(J1(:));
%                 J1 = J1 / (max(J1(:)) - min(J1(:)));
%                 
%                 J2 = mean(Seg,3);
%                 J2 = J2 - min(J2(:));
%                 J2 = J2 / (max(J2(:)) - min(J2(:)));
  
                J1 = mean(Jrecon,3);
                clim = quantile(J1(:),qlim);
                J1 = (J1 - clim(1))/diff(clim);
                
                J2 = mean(Seg,3);
                clim = quantile(J2(:),qlim);
                J2 = (J2 - clim(1))/diff(clim);
                
                COV = cov(J1(:),J2(:));
                J2_ = (J2 - mean(J2(:)))/COV(2,2)*COV(1,2) + mean(J1(:));
                
                imagesc(xJ{f},yJ{f},cat(3,J1,J2_,J1))
                
            end
            title(['z = ' num2str(zJ(f))])
            axis image;
            xlabel x;
            ylabel y;
            set(gca,'xlim',6000*[-1,1],'ylim',5000*[-1 1]);
            if  ~isempty(strfind(thename_,'rat'))
                % update for rat
                
                set(gca,'xlim',9000*[-1,1],'ylim',8000*[-1 1]);

            end
            [directory_, filename_, extension_] = fileparts(files{f});
            saveas(5,[outdir filename_ '_preview_' num2str(t) '_straight.png'])
            
            
            
            
            %% 
            % make a geojson for level lines            
            % we will sample at 1mm
            % and use 50micron resolution
            % three files, for phiiAxyziPhiJiAJiY, phiiAxyziPhiJiAJiX, phiiAxyziPhiJiAJiZ
            levels = -6000:500:6000;
            contoursx = {};
            contoursy = {};
            contoursz = {};
            for l = 1 : length(levels)
            contoursx{l} = contourc(xJ{f},yJ{f},phiiAxyziPhiJiAJiX,levels(l)*[1,1]);
            contoursy{l} = contourc(xJ{f},yJ{f},phiiAxyziPhiJiAJiY,levels(l)*[1,1]);
            contoursz{l} = contourc(xJ{f},yJ{f},phiiAxyziPhiJiAJiZ,levels(l)*[1,1]);
            end
            clf
            imagesc(xJ{f},yJ{f},Jrecon)
            hold on;
            contour(xJ{f},yJ{f},phiiAxyziPhiJiAJiX,levels,'r');
            contour(xJ{f},yJ{f},phiiAxyziPhiJiAJiY,levels,'g');
            contour(xJ{f},yJ{f},phiiAxyziPhiJiAJiZ,levels,'b');
            hold off
            % do x contours
            contoursx_ = {};            
            namesx_ = {};
            contoursy_ = {};            
            namesy_ = {};
            contoursz_ = {};            
            namesz_ = {};
            for l = 1 : length(levels)
                namesx_{l} = num2str(levels(l));
                sc = contoursx{l};   
                count = 0;
                while 1
                    if isempty(sc)
                        break
                    end
                    count = count + 1;
                    n = sc(2,1);

                    xdata = sc(1,2:n+1);
                    ydata = sc(2,2:n+1);

                    contoursx_{l}{count} = [xdata;ydata];                    
                    sc = sc(:,n+2:end);

                end
            end
            empty = cellfun(@(x) isempty(x),contoursx_);
            contoursx_ = contoursx_(~empty);
            namesx_ = namesx_(~empty);
            
            for l = 1 : length(levels)
                namesy_{l} = num2str(levels(l));
                sc = contoursy{l};                
                count = 0;
                while 1
                    if isempty(sc)
                        break
                    end
                    count = count + 1;
                    n = sc(2,1);

                    xdata = sc(1,2:n+1);
                    ydata = sc(2,2:n+1);

                    contoursy_{l}{count} = [xdata;ydata];                    
                    sc = sc(:,n+2:end);

                end
            end
            empty = cellfun(@(x) isempty(x),contoursy_);
            contoursy_ = contoursy_(~empty);
            namesy_ = namesy_(~empty);
            for l = 1 : length(levels)
                namesz_{l} = num2str(levels(l));
                sc = contoursz{l};  
                count = 0;
                while 1
                    if isempty(sc)
                        break
                    end
                    count = count + 1;
                    n = sc(2,1);

                    xdata = sc(1,2:n+1);
                    ydata = sc(2,2:n+1);

                    contoursz_{l}{count} = [xdata;ydata];                    
                    sc = sc(:,n+2:end);

                end
            end
            empty = cellfun(@(x) isempty(x),contoursz_);
            contoursz_ = contoursz_(~empty);
            namesz_ = namesz_(~empty);

            % write out the json
            contournamex = [outdir filename_  '_grid_x_straight.geojson'];
            contournamey = [outdir filename_  '_grid_y_straight.geojson'];
            contournamez = [outdir filename_  '_grid_z_straight.geojson'];
            
            write_geojson_grid(contournamex,namesx_,contoursx_);
            write_geojson_grid(contournamey,namesy_,contoursy_);
            write_geojson_grid(contournamez,namesz_,contoursz_);
            
            
            
        end
    end
    
    
    %% 
    % output the transformation that was used to deform the labels
    % this is registered to atlas
    out = single(cat(4,phiiAxyziPhiJiAJiX-XJ,phiiAxyziPhiJiAJiY-YJ, phiiAxyziPhiJiAJiZ-zJ(f)));
    write_vtk_image(xJ{f},yJ{f},zJ_write,out,[outdir 'registered_to_atlas_displacement_' fname_ '.vtk'],'registered_to_atlas_displacement');
    
    
    
    
    %%
    % last is we need a grid
    % I will use isocontours on phiiAxyziPhiJiAJiX
    % note Y and Z are important
%     hold on;
%     contour(xJ{f},yJ{f},phiiAxyziY,[-6000 : 500 : 6000],'k');
%     contour(xJ{f},yJ{f},phiiAxyziZ,[-6000 : 500 : 6000],'k');
%     contour(xJ{f},yJ{f},phiiAxyziX,[-6000 : 500 : 6000],'k'); 
% % note deformation out of plane spans many slices
%     hold off;
    


    
    %% 
    % now do it again with the detjac instead of the seg
    F = griddedInterpolant({yV,xV,zV},detjac,'linear','nearest');
    distortion = F(phiiAxyziPhiJiAJiY, phiiAxyziPhiJiAJiX, phiiAxyziPhiJiAJiZ);
    write_vtk_image(xJ{f},yJ{f},zJ,single(distortion),[outdir 'atlas_to_registered_detjac_' fname_ '.vtk'],'atlas_to_registered_detjac')

    
    %%
    drawnow;

    
end % of loop over files



%%
% structure the output directory
keyboard