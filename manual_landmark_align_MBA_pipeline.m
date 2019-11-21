function manual_landmark_align_MBA_pipeline(input_dir, qc_file)
% model can be
% rigid
% affine
% spline
% rigid+spline
% affine+spline


if nargin == 0
    input_dir = '/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/Registration/PMD1027_case1_data/INPUT_DATA/';
    qc_file = '/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/Registration/PMD1027_case1_data/OUTPUT/PMD1028&1027-F36-2013.03.02-02.15.31_PMD1027_3_0108_preview_1_straight.png';
end


% for Xu, the "detailed output" directory is the same as the output
% directory
% detailed_postfix = '_detailed';
detailed_postfix = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parse geometry file
geometry_file = dir([input_dir '*.csv']);
fid = fopen([input_dir geometry_file(1).name],'rt');
line = fgetl(fid); % ignore the first line
% filename, nx, ny, nz, dx, dy, dz, x0, y0, z0
csv_data = {};
nissl_pattern = '*-N*.tif';
fluoro_pattern = '*-F*.tif';
count = 0;
while 1
    line = fgetl(fid);
    if line == -1
        break
    end
    % check pattern
    count = count + 1;
    
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






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get filename from qc
[dir_,fname_,ext_] = fileparts(qc_file);
ind = strfind(fname_,'_preview');
num = fname_(ind-4:ind-1);
inds = ~cellfun(@(x)isempty(x),strfind(files,['_' num '.tif']));


inds = zeros(1,length(zJ0));
for i = 1 : length(zJ0)
    if ~is_fluoro(i)
        continue
    end
    
    if isempty(strfind(files{i},['_' num '.tif']))
        continue
    end
    
    fluoro_ind = i;
    
    
    
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
    
    
    
    
end





addpath Functions/plotting
figure(111);

% we will add landmarks in a loop
% these are the choices
% add to left
% add to right
% undo
% exit
QI = zeros(0,3);
QJ = zeros(0,3);

A = eye(3);
[XJ,YJ] = meshgrid(xJ,yJ);

while 1
    
    hI = subplot(2,2,1);
    imagesc(xI,yI,I);   
    hold on;
    scatter(QI(:,1),QI(:,2));
    hold off;
    axis image;
    
    
    hJ = subplot(2,2,2);
    imagesc(xJ,yJ,J);
    hold on;
    scatter(QJ(:,1),QJ(:,2));
    hold off;
    
    axis image
    
    
    title(hI,'Select point')
    title(hJ,'')
    [qxI,qyI,BUTTON] = ginput(1);
    if isempty(qxI) % if you hit enter, its over
        break
    end

    if BUTTON == 2
        break
        continue
    end
    if BUTTON == 3
        if size(QI,1) >= 1
            QI = QI(1:end-1,:);
        end
        QJ = QJ(1:size(QI,1),:);
        continue
    end
    
    
    title(hJ,'Select point')
    title(hI,'')

    [qxJ,qyJ,BUTTON] = ginput(1);
    if isempty(qxJ)
        break
    end
    if BUTTON == 2
        break
        continue
    end
    if BUTTON == 3
        if size(QI,1) >= 1
            QI = QI(1:end-1,:);
        end
        QJ = QJ(1:size(QI,1),:);
        continue
    end
    % update tform
    QI = [QI;qxI,qyI,1];
    QJ = [QJ;qxJ,qyJ,1];
    
    if size(QI,1) < 3
        model = 'translation';
    else
        model = 'rigid';        
    end
    if strfind(model,'translation')
        A = eye(3);
        A(1:2,end) = mean(QJ(:,1:2),1) - mean(QI(:,1:2),1);
    end
    if strfind(model,'affine')
        A = ((QI'*QI)\(QI'*QJ))';
    end
    if strfind(model,'rigid')
        T = mean(QJ,1) - mean(QI,1);
        QI0 = bsxfun(@minus,QI,mean(QI,1));
        QJ0 = bsxfun(@minus,QJ,mean(QJ,1));
        M = QI0(:,1:2)' * QJ0(:,1:2);
        [U,S,V] = svd(M);
        % V' R' U should be equal to S
        % we want to maximize
        % trace(RX Yt)
        % = trace(R M) for M = X Yt
        % = trace(R USV')
        % = trace ([V' R U] S)
        % that means
        % V' R U should be diagonal
        % and it should have diagonal values that line up with the sign of
        % S
        R = U * sign(S) * V';
        A = [R',T(1:2)';0,0,1];
        
        
    end
    
    AI = zeros(size(J));
    B = inv(A);
    Xs = B(1,1)*XJ + B(1,2)*YJ + B(1,3);
    Ys = B(2,1)*XJ + B(2,2)*YJ + B(2,3);

    
    
    if strfind(model,'spline')
        % now spline matching on the residuals
        % this is gonna be harder and slower
        % to deform I to match J, I need to match J points to I points
        BQJ = QJ * B';
        Delta = QI - BQJ;
        width = 100;
        noise = 10;
        K = exp(-pdist2(BQJ,BQJ).^2 / 2 / width.^2);
        P = (K + eye(size(QI,1))*noise^2)\Delta;
        VX = zeros(size(J,1),size(J,2));
        VY = zeros(size(J,1),size(J,2));
        for i = 1 : size(QI,1)
            K = exp( - ((XJ - BQJ(i,1)).^2 + (YJ - BQJ(i,2)).^2)/2/width.^2 );
            VX = VX + P(i,1)*K;
            VY = VY + P(i,2)*K;
        end
        Xs = Xs + VX;
        Ys = Ys + VY;
        
    end
    
    
    for c = 1 : size(I,3)
        F = griddedInterpolant({yI,xI},I(:,:,c),'linear','nearest');
        AI(:,:,c) = F(Ys,Xs);
    end
    
    subplot(2,2,3)
    imagesc(xJ,yJ,AI*0.5 + J*0.5);
    axis image
    title(['enter->finish, middl click->auto reg'])
    
    
    
    drawnow
    
    
end



% final update (just copied from above
if strfind(model,'translation')
    A = eye(3);
    A(1:2,end) = mean(QJ(:,1:2),1) - mean(QI(:,1:2),1);
end

if strfind(model,'affine')
    A = ((QI'*QI)\(QI'*QJ))';
end
if strfind(model,'rigid')
    T = mean(QJ,1) - mean(QI,1);
    QI0 = bsxfun(@minus,QI,mean(QI,1));
    QJ0 = bsxfun(@minus,QJ,mean(QJ,1));
    M = QI0(:,1:2)' * QJ0(:,1:2);
    [U,S,V] = svd(M);
    % V' R' U should be equal to S
    % we want to maximize
    % trace(RX Yt)
    % = trace(R M) for M = X Yt
    % = trace(R USV')
    % = trace ([V' R U] S)
    % that means
    % V' R U should be diagonal
    % and it should have diagonal values that line up with the sign of
    % S
    R = U * sign(S) * V';
    A = [R',T(1:2)';0,0,1];
    
    
end

AI = zeros(size(J));
B = inv(A);
Xs = B(1,1)*XJ + B(1,2)*YJ + B(1,3);
Ys = B(2,1)*XJ + B(2,2)*YJ + B(2,3);



if strfind(model,'spline')
    % now spline matching on the residuals
    % this is gonna be harder and slower
    % to deform I to match J, I need to match J points to I points
    BQJ = QJ * B';
    Delta = QI - BQJ;
    width = 100;
    noise = 10;
    K = exp(-pdist2(BQJ,BQJ).^2 / 2 / width.^2);
    P = (K + eye(size(QI,1))*noise^2)\Delta;
    VX = zeros(size(J,1),size(J,2));
    VY = zeros(size(J,1),size(J,2));
    for i = 1 : size(QI,1)
        K = exp( - ((XJ - BQJ(i,1)).^2 + (YJ - BQJ(i,2)).^2)/2/width.^2 );
        VX = VX + P(i,1)*K;
        VY = VY + P(i,2)*K;
    end
    Xs = Xs + VX;
    Ys = Ys + VY;
    
end


for c = 1 : size(I,3)
    F = griddedInterpolant({yI,xI},I(:,:,c),'linear','nearest');
    AI(:,:,c) = F(Ys,Xs);
end

subplot(2,2,3)
imagesc(xJ,yJ,AI*0.5 + J*0.5);
axis image

drawnow
%%
% now run the automatic
if BUTTON == 2
    downs = [16,8,6,3]; % note that 3 takes us from 15 to 45, the res we're working at
    niter = 100;
    NtoF = slice_to_slice_rigid_alignment_GN_weight(xI,yI,I,xJ,yJ,J,A,downs,niter);
else
    NtoF = A;
end

%%
% now insert into pipeline
NtoFfile = [dir_ detailed_postfix '/NtoF.mat'];
varsNtoF = load(NtoFfile);
backup = [NtoFfile(1:end-4) '_beforeQC.mat'];
if ~exist(backup,'file')
    copyfile(NtoFfile,backup);
end


matfile = [dir_ detailed_postfix '/combined_A.mat'];

% make a copy if it doesn't exist
backup = [matfile(1:end-4) '_beforeQC.mat'];
if ~exist(backup,'file')
    copyfile(matfile,backup)
end
vars = load(matfile);

% restore the original by inverting the "bad" NtoF
AJorig = varsNtoF.NtoF(:,:,fluoro_ind) \ vars.AJ(:,:,fluoro_ind);

% update and save
AJ = vars.AJ;
A = vars.A;
AJ(:,:,fluoro_ind) = NtoF * AJorig;
save(matfile,'A','AJ');

% save the new N to F (in case I need to invert it again)
varsNtoF.NtoF(:,:,fluoro_ind) = NtoF;
is_nissl = varsNtoF.is_nissl;
is_fluoro = varsNtoF.is_fluoro;
files = varsNtoF.files;
NtoF = varsNtoF.NtoF;
inds = varsNtoF.inds;
save(NtoFfile,'is_nissl','is_fluoro','files','NtoF','inds')

%%
% after you are finished all the slices, run this code
% seg_file = '/cis/home/dtward/Documents/ARA/Mouse_CCF/vtk/annotation_50.vtk';
% atlas_file = '/cis/home/dtward/Documents/ARA/Mouse_CCF/vtk/ara_nissl_50.vtk';
% apply_deformation({seg_file,atlas_file}, input_dir, [dir_ detailed_postfix '/'], dir_);

