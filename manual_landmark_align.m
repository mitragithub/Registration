function [Xs,Ys] = manual_landmark_align(I,J,model)
% model can be
% rigid
% affine
% spline
% rigid+spline
% affine+spline


if nargin == 0
    J = imread(['/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/Xu2Daniel/PMD3317/PMD3317&3316-F64-2019.02.25-22.56.20_PMD3317_3_0192.tif']);
    I = imread(['/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/Xu2Daniel/PMD3317/PMD3317&3316-N64-2019.02.06-04.25.48_PMD3317_3_0192.tif']);
    
    I = double(I)/255.0;
    J = double(J)/255.0;
    model = 'rigid';
%     model = 'affine';
%     model = 'rigid+spline';
%     model = 'affine+spline';
%     keyboard
end


figure;

% we will add landmarks in a loop
% these are the choices
% add to left
% add to right
% undo
% exit
QI = zeros(0,3);
QJ = zeros(0,3);

A = eye(3);
xI = 1 : size(I,2); xI = xI - mean(xI);
yI = 1 : size(I,1); yI = yI - mean(yI);

xJ = 1 : size(J,2); xJ = xJ - mean(xJ);
yJ = 1 : size(J,1); yJ = yJ - mean(yJ);

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
    
    
end


% final update (just copied from above
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