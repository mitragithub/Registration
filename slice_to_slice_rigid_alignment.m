function A = slice_to_slice_rigid_alignment(xI,yI,I,xJ,yJ,J,A0,downs,niter,eTfactor,eLfactor)
% rigidly align two slices, possibly different modalities
% we wish to align image I to image J rigidly
% we use exponential parameterization
% minimize
% int 1/2 |f(I(Ai x)) - J(x))|^2 dx
% let fI = I'
% int 1/2 |I'(Ai x) - J(x))|^2 dx
% right perturbation
% A \to A exp(e dA)
% or Ai \to  exp(-edA)Ai
% d_ee int 1/2 |I'(exp(-edA) Ai x) - J(x))|^2 dx |_e=0
% = \int (I'(Ai x) - J(x))^T DI'(Ai x)  (-1) dA Ai x dx
%= (-1) \int (I'(Ai x) - J(x))^T D[I'(Ai x)] A  dA Ai x dx
if nargin == 0
    % run an example
    
    I = imread('/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/Xu2Daniel/PMD964/PMD965&964-N33-2012.06.21-00.30.07_PMD964_3_0099.tif');
    I = double(I)/255.0;
    dxI = [14.72 14.72];
    xI = (0:size(I,2)-1)*dxI(1);xI = xI - mean(xI);
    yI = (0:size(I,1)-1)*dxI(2);yI = yI - mean(yI);
    
    J = imread('/cis/home/dtward/Documents/intensity_transform_and_missing_data/csh_slices/Xu2Daniel/PMD964/PMD965&964-F33-2012.06.17-14.22.42_PMD964_3_0099.tif');
    J = double(J)/255.0;
    dxJ = [14.72 14.72];
    xJ = (0:size(J,2)-1)*dxJ(1);xJ = xJ - mean(xJ);
    yJ = (0:size(J,1)-1)*dxJ(2);yJ = yJ - mean(yJ);
    
    order = 2;
    keyboard
end
addpath Functions/plotting
addpath Functions/gradient
addpath Functions/downsample

% I = mean(I,3);
% J = mean(J,3);

danfigure(1);
imagesc(xI,yI,I)
axis image

danfigure(2);
imagesc(xJ,yJ,J)
axis image

CI = size(I,3);
CJ = size(J,3);

[XJ,YJ] = meshgrid(xJ,yJ);
dxI = [xI(2) - xI(1), yI(2) - yI(1)];
dxJ = [xJ(2) - xJ(1), yJ(2) - yJ(1)];

if nargin < 7
    A0 = eye(3);
end
if nargin < 8
    downs = [8,4,2];
end
if nargin < 9
    niter = 100;
end
if length(niter) == 1
    niter = repmat(niter,length(downs));
end
% step sizes to multiply by d
if nargin < 10
    eTfactor = 1e-2;
end
if nargin < 11
    eLfactor = 1e-8;
end
A = A0;
% downsample

ndraw = 11;
for downloop = 1 : length(downs)
    d = downs(downloop);
    eL = eLfactor*d;
    eT = eTfactor*d;
    e = [1,1,0;
        1,1,0;
        0,0,0]*eL + ...
        [0,0,1;
        0,0,1;
        0,0,0]*eT;
    Id = [];
    for c = 1 : CI
        [~,~,Id(:,:,c)] = downsample2D(1:size(I,2),1:size(I,1),I(:,:,c),[1,1]*d);
    end
    dxId = dxI*d;
    xId = (1 : size(Id,2))*dxId(1); xId = xId - mean(xId);
    yId = (1 : size(Id,1))*dxId(2); yId = yId - mean(yId);
    Jd = [];
    for c = 1 : CI
        [~,~,Jd(:,:,c)] = downsample2D(1:size(J,2),1:size(J,1),J(:,:,c),[1,1]*d);
    end
    dxJd = dxJ*d;
    xJd = (1 : size(Jd,2))*dxJd(1); xJd = xJd - mean(xJd);
    yJd = (1 : size(Jd,1))*dxJd(2); yJd = yJd - mean(yJd);
    
    
    
    %
    [XJd,YJd] = meshgrid(xJd,yJd);
    for it = 1 : niter(downloop)
        
        % transform I
        Ai = inv(A);
        AId = zeros(size(Jd,1),size(Jd,2),CI);
        Xs = Ai(1,1)*XJd + Ai(1,2)*YJd + Ai(1,3);
        Ys = Ai(2,1)*XJd + Ai(2,2)*YJd + Ai(2,3);
        for c = 1 : CI
            F = griddedInterpolant({yId,xId},Id(:,:,c),'linear','nearest');
            AId(:,:,c) = F(Ys,Xs);
        end
        
        
        
        % now intensity
        % start with just linear
        D = cat(3,ones(size(AId,1),size(AId,2)),AId);
        Dvec = reshape(D,[],size(D,3));
        if it == 1
            coeffs = (Dvec'*Dvec) \ (Dvec' * reshape(Jd,[],CJ));
        end
        fAId = reshape(Dvec*coeffs,size(Jd));
        
        
        
        % now get the cost
        err = fAId - Jd;
        
        
        if ~mod(it-1,ndraw) || it == niter(downloop)
        danfigure(1);
        imagesc(xId,yId,Id)
        axis image
        title('I')
        
        danfigure(2);
        imagesc(xJd,yJd,Jd)
        axis image
        title('J')
        
        danfigure(3);
        imagesc(xJd,yJd,AId);
        axis image
        title('AI')
        
        danfigure(4);
        imagesc(xJd,yJd,fAId);
        axis image
        title('fAI')
        danfigure(5);
        imagesc(xJd,yJd,err/2.0 + 0.5)
        axis image
        title('err')
        end
        
        E = sum(err(:).^2)/2*prod(dxJd);
        disp(['Iteration ' num2str(it) '/' num2str(niter(downloop)) ', energy ' num2str(E)]);
        grad = zeros(3,3);
        % and the gradient
        fAId_x = zeros(size(fAId));
        fAId_y = zeros(size(fAId));
        for c = 1 : CJ
            [fAId_x(:,:,c),fAId_y(:,:,c)] = gradient2d(fAId(:,:,c),dxJd(1),dxJd(2));
        end
        for r = 1 : 2
            for c = 1 : 3
                dA = double((1:3==r))' * double((1:3==c));
                AdAAi = A*dA/A;
                Xs = AdAAi(1,1)*XJd + AdAAi(1,2)*YJd + AdAAi(1,3);
                Ys = AdAAi(2,1)*XJd + AdAAi(2,2)*YJd + AdAAi(2,3);
                integrand = err.*(bsxfun(@times, fAId_x,Xs) + bsxfun(@times, fAId_y,Ys));
                grad(r,c) = -sum(integrand(:)).*prod(dxJd);
            end
        end
        % rigid
        grad(1:2,1:2) = grad(1:2,1:2) - grad(1:2,1:2)';
        
        % update
        A = A * expm(-e.*grad);
        
        
        
        
        drawnow
        
        
    end % of iter
    
end % of downsample