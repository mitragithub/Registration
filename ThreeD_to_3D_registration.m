function [phiiAix,phiiAiy,phiiAiz,Aphix,Aphiy,Aphiz,A,vtx,vty,vtz] = ThreeD_to_3D_registration(...
    xI,yI,zI,I,...
    xJ,yJ,zJ,J,...
    A0,...
    nt,a,p,...
    sigmaR,sigmaM,sigmaA,sigmaB,prior,...
    order,...
    eA,eV,apfactor,...
    downs,niter,...
    OPT)
%
% joint estimation of affine and diffeomorphic registration between 3D
% datasets
% TODO: make IO consistent with other code


addpath Functions/downsample
addpath Functions/plotting
addpath Functions/gradient
if nargin == 0
    addpath Functions/vtk
    
    % get an example dataset for debugging
    [xI,yI,zI,I,title_,names] = read_vtk_image('MD787_atlas_vtk/CSHL_A_ExVivo_37_SE_TE008ms_50micron.vtk');
    [xJ,yJ,zJ,J,title_,names] = read_vtk_image('atlas_50_vtk/ara_nissl_50.vtk');
    
    % squash bright signals
    % what fraction of the voxels are 0
    frac = sum(J(:)>0)/numel(J); % 0.4877, almost half
    % I want to normalize the ones that are greater than 0
    J(J>0) = tiedrank(J(J>0));
    
    
    I(I>0.2) = tiedrank(I(I>0.2));
    
    A0 = eye(4);

    downs = [4,2];
    niter = 20;
    eV = 1e5; % velocity field step size
    eA = 0.5; %
    
    
    % other params
    nt = 5; % number of timesteps
    a = 100; % 100 microns
    apfactor = 2; % 2 voxels
    p = 2; % power of linaer operators
    
    sigmaR = 1e6;
    sigmaM = 1.0;
    
    % contrast
    order = 3;
    
    % weights
    prior = [0.8,0.1,0.1];
    sigmaA = 5.0;
    sigmaB = 2.0;
    
    
end

if nargin < 24
    OPT = struct;
end

% EXPERIMENTAL gauss newton approximation feature for velocity
velocity_GN = 0;
if isfield(OPT,'velocity_GN')
    velocity_GN = OPT.velocity_GN;
end

if length(niter) == 1
    niter = ones(size(downs))*niter;
end


I = (I - mean(I(:))) / std(I(:));
J = (J - mean(J(:))) / std(J(:));

% get ready for downsampling
I0 = I;
xI0 = xI;
yI0 = yI;
zI0 = zI;

J0 = J;
xJ0 = xJ;
yJ0 = yJ;
zJ0 = zJ;


% for plotting
nplot = 5;
qlim = [0.001,0.999];
climI = quantile(I(:),qlim);
climJ = quantile(J(:),qlim);

% initialize
A = A0;
Ai = inv(A);

% contrasts
muA = max(J(:));
muB = min(J(:));
PI = ones(3,1)/3;
    
    
    
E = zeros(1,sum(niter));
% loop over downsamples
itercount = 0;
for downloop = 1 : length(downs)
    down = downs(downloop);
    
    % downsample the images
    [xI,yI,zI,I] = downsample(xI0,yI0,zI0,I0,[1,1,1]*down);
    [XI,YI,ZI] = meshgrid(xI,yI,zI);
    dxI = [xI(2)-xI(1), yI(2)-yI(1), zI(2)-zI(1)];
    nxI = [size(I,2),size(I,1),size(I,3)]; 
    fxI = (0:nxI(1)-1)/nxI(1)/dxI(1);
    fyI = (0:nxI(2)-1)/nxI(2)/dxI(2);
    fzI = (0:nxI(3)-1)/nxI(3)/dxI(3);
    [FXI,FYI,FZI] = meshgrid(fxI,fyI,fzI);
    
    
    danfigure(1);
    sliceView(xI,yI,zI,I,nplot,climI);
    title('atlas')
    
    [xJ,yJ,zJ,J] = downsample(xJ0,yJ0,zJ0,J0,[1,1,1]*down);
    dxJ = [xJ(2)-xJ(1), yJ(2)-yJ(1), zJ(2)-zJ(1)];        
    [XJ,YJ,ZJ] = meshgrid(xJ,yJ,zJ);
    % we need a padded version
    xJp = [xJ(1)-dxJ(1),xJ,xJ(end)+dxJ(1)];
    yJp = [yJ(1)-dxJ(2),yJ,yJ(end)+dxJ(2)];
    zJp = [zJ(1)-dxJ(3),zJ,zJ(end)+dxJ(3)];
    danfigure(2);
    sliceView(xJ,yJ,zJ,J,nplot,climJ);
    title('target')
    
    

    
    % initialze velocity at this res
    if downloop == 1
        vtx = zeros(size(I,1),size(I,2),size(I,3),nt);
        vty = zeros(size(I,1),size(I,2),size(I,3),nt);
        vtz = zeros(size(I,1),size(I,2),size(I,3),nt);
    else
        vtx_ = vtx;
        vty_ = vty;
        vtz_ = vtz;
        vtx = zeros(size(I,1),size(I,2),size(I,3),nt);
        vty = zeros(size(I,1),size(I,2),size(I,3),nt);
        vtz = zeros(size(I,1),size(I,2),size(I,3),nt);
        for t = 1 : nt
            vtx(:,:,:,t) = upsample(vtx_(:,:,:,t),size(I));
            vty(:,:,:,t) = upsample(vty_(:,:,:,t),size(I));
            vtz(:,:,:,t) = upsample(vtz_(:,:,:,t),size(I));
        end
    end

    % initialze weights at this res
    if downloop == 1
        WM = ones(size(J))*prior(1);
        WA = ones(size(J))*prior(2);
        WB = ones(size(J))*prior(3);
    else        
        WM = upsample(WM,size(J));
        WA = upsample(WA,size(J));
        WB = upsample(WB,size(J));        
    end

    
    % create kernels at this res
    LL = ( 1.0 - a^2*2*( (cos(2.0*pi*FXI*dxI(1))-1.0)/dxI(1)^2 + (cos(2.0*pi*FYI*dxI(2))-1.0)/dxI(2)^2 + (cos(2.0*pi*FZI*dxI(3))-1.0)/dxI(3)^2  ) ).^(2*p);
    K = 1.0./LL;
    % Kp for preconditioner
    ap = sqrt(sum(dxI.^2))*apfactor;
    LLp = ( 1.0 - (ap)^2*2*( (cos(2.0*pi*FXI*dxI(1))-1.0)/dxI(1)^2 + (cos(2.0*pi*FYI*dxI(2))-1.0)/dxI(2)^2 + (cos(2.0*pi*FZI*dxI(3))-1.0)/dxI(3)^2  ) ).^(2*p);
    Kp = LL./LLp;

    % vmax should be half a voxel
    Vmax = sqrt(sum(dxI.^2))*0.5*100;
    
    
    
    % start gradient descent loop    
    for it = 1 : niter(downloop)
        itercount = itercount + 1;
        
        % build diffeo
        phiix = XI;
        phiiy = YI;
        phiiz = ZI;
        It = zeros(size(I,1),size(I,2),size(I,3),nt);
        It(:,:,:,1) = I;
        dt = 1/nt;
        for t = 1 : nt
            % deform the image
            if t > 1 
                F = griddedInterpolant({yI,xI,zI},I,'linear','nearest');
                It(:,:,:,t) = F(phiiy,phiix,phiiz);
            end
            % update the phii
            vx = vtx(:,:,:,t);vy = vty(:,:,:,t);vz = vtz(:,:,:,t);
            Xs = XI - dt*vx; Ys = YI - dt*vy; Zs = ZI - dt*vz;
            F = griddedInterpolant({yI,xI,zI},phiix-XI,'linear','nearest');
            phiix = F(Ys,Xs,Zs) + Xs;
            F = griddedInterpolant({yI,xI,zI},phiiy-YI,'linear','nearest');
            phiiy = F(Ys,Xs,Zs) + Ys;
            F = griddedInterpolant({yI,xI,zI},phiiz-ZI,'linear','nearest');
            phiiz = F(Ys,Xs,Zs) + Zs;
                        
        end
        
        % cost
        vtxhat = fft(fft(fft(vtx,[],1),[],2),[],3);
        vtyhat = fft(fft(fft(vty,[],1),[],2),[],3);
        vtzhat = fft(fft(fft(vtz,[],1),[],2),[],3);
        
        % combine with affine
        
        
        % transform I
        Xs = Ai(1,1)*XJ + Ai(1,2)*YJ + Ai(1,3)*ZJ + Ai(1,4);
        Ys = Ai(2,1)*XJ + Ai(2,2)*YJ + Ai(2,3)*ZJ + Ai(2,4);
        Zs = Ai(3,1)*XJ + Ai(3,2)*YJ + Ai(3,3)*ZJ + Ai(3,4);
        
        F = griddedInterpolant({yI,xI,zI},phiix-XI,'linear','nearest');
        phiiAix = F(Ys,Xs,Zs) + Xs;
        F = griddedInterpolant({yI,xI,zI},phiiy-YI,'linear','nearest');
        phiiAiy = F(Ys,Xs,Zs) + Ys;
        F = griddedInterpolant({yI,xI,zI},phiiz-ZI,'linear','nearest');
        phiiAiz = F(Ys,Xs,Zs) + Zs;
        
        F = griddedInterpolant({yI,xI,zI},I,'linear','nearest');
        AphiI = F(phiiAiy,phiiAix,phiiAiz);
        
        
        danfigure(3);
        sliceView(xJ,yJ,zJ,AphiI,nplot,climI)
        
%         % linear prediction
%         COV = cov(AphiI(:),J(:));
%         fAI = (AphiI - mean(AphiI(:)))/COV(1,1)*COV(1,2) + mean(J(:));
        
        % polynomial prediction
        % basis
        B = ones(numel(J),order + 1);
        DB = zeros(numel(J),order + 1); DB(:,2) = 1;
        for o = 2 : order + 1
            B(:,o) = AphiI(:).^(o-1);
            if o > 2
                DB(:,o) = (o-1) * AphiI(:).^(o-2);
            end
        end
        % find coefficients
        coeffs = (B'*B)\(B' * J(:));
        fAphiI = reshape(B * coeffs,size(J));
        danfigure(4);
        sliceView(xJ,yJ,zJ,fAphiI,nplot,climJ)
        
        % error
        err = (fAphiI - J);
        danfigure(5);
        sliceView(xJ,yJ,zJ,err,nplot,[-1,1]*max(abs(climJ)));
        
        % cost
        EM(itercount) = sum(abs(err(:)).^2.*WM(:))/2*prod(dxJ);
        ER(itercount) = sum(sum(sum(sum( bsxfun(@times, (abs(vtxhat).^2 + abs(vtyhat).^2 + abs(vtzhat).^2) ,  LL ) ))))*prod(dxI)*dt/numel(vtx(:,:,:,1))/sigmaR^2/2;
        E(itercount) = EM(itercount) + ER(itercount);
        danfigure(6);
        plot(ER(1:itercount))
        hold on;
        plot(EM(1:itercount))
        plot(E(1:itercount))
        hold off;
        disp(['Down loop ' num2str(downloop) ', iteration '  num2str(it) ', energy ' num2str(E(itercount)) ' (match ' num2str(EM(itercount)) ', reg ' num2str(ER(itercount)) ')'])
        
        
        
        % update weights
        if ~mod(itercount-1,2)
%             PI = [1,1,1]/3; % not using this weighting
            WM = 1/sqrt(2.0*pi*sigmaM^2)*exp(-(fAphiI-J).^2/2.0/sigmaM^2)*prior(1)*PI(1);
            WA = 1/sqrt(2.0*pi*sigmaA^2)*exp(-(muA-J).^2/2.0/sigmaA^2)*prior(2)*PI(2);
            WB = 1/sqrt(2.0*pi*sigmaB^2)*exp(-(muB-J).^2/2.0/sigmaB^2)*prior(3)*PI(3);
            WSum = WM+WA+WB;
            WSum = WSum + max(WSum(:))*1e-6;
            WM = WM./WSum;
            WA = WA./WSum;
            WB = WB./WSum;
            danfigure(7);
            sliceView(xJ,yJ,zJ,cat(4,WM,WA,WB),nplot,[0,1]);

        end
        
        % now we take gradient
        % to do this use right perturbation
        % int (I(Ai x) - J(x))^2 dx
        % d_de int (I(exp(-e dA )Ai x) - J(x))^2 dx e = 0
        % = int (I(Ai x) - J(x)) DI (Ai x) (-1) dA Ai x dx
        % = int (I(Ai x) - J(x)) D[I (Ai x)] A (-1) dA Ai x dx
        
        [fAphiI_x,fAphiI_y,fAphiI_z] = gradient3d(fAphiI,dxJ(1),dxJ(2),dxJ(3));
        
        Jerr = zeros(size(J,1),size(J,2),size(J,3),12);
        count = 0;
        for r = 1 : 3
            for c = 1 : 4
                dA = double((1:4==r))' * double((1:4==c));
                AdAAi = A*dA;
                Xs = AdAAi(1,1)*XJ + AdAAi(1,2)*YJ + AdAAi(1,3)*ZJ + AdAAi(1,4);
                Ys = AdAAi(2,1)*XJ + AdAAi(2,2)*YJ + AdAAi(2,3)*ZJ + AdAAi(2,4);
                Zs = AdAAi(3,1)*XJ + AdAAi(3,2)*YJ + AdAAi(3,3)*ZJ + AdAAi(3,4);
                count = count + 1;
                Jerr(:,:,:,count) = (bsxfun(@times, fAphiI_x,Xs) + bsxfun(@times, fAphiI_y,Ys) + bsxfun(@times, fAphiI_z,Zs));
            end
        end
        JerrJerr = squeeze(sum(sum(sum(bsxfun(@times, permute(bsxfun(@times,Jerr,WM),[1,2,3,4,5]) , permute(Jerr,[1,2,3,5,4])),3),2),1));
        
        % step
        Aistep = JerrJerr \ squeeze(sum(sum(sum(bsxfun(@times, Jerr, bsxfun(@times,err,WM)),3),2),1));
        Aistep = reshape(Aistep,4,3)';
        
        
        
        
        
        % get the error to pull back
        errDf = err.*WM .* reshape(DB*coeffs,size(err));
        errDfp = zeros(size(errDf)+2);
        errDfp(2:end-1,2:end-1,2:end-1) = errDf;
        if velocity_GN
            WMp = zeros(size(errDf)+2);
            WMp(2:end-1,2:end-1,2:end-1) = WM;
        end
        
        % build diffeo backwards
        phix = XI;
        phiy = YI;
        phiz = ZI;
        
        dt = 1/nt;
        for t = nt : -1 : 1
            % deform the image
            
            % update the phii
            vx = vtx(:,:,:,t);vy = vty(:,:,:,t);vz = vtz(:,:,:,t);
            Xs = XI + dt*vx; Ys = YI + dt*vy; Zs = ZI + dt*vz; %note + sign
            F = griddedInterpolant({yI,xI,zI},phix-XI,'linear','nearest');
            phix = F(Ys,Xs,Zs) + Xs;
            F = griddedInterpolant({yI,xI,zI},phiy-YI,'linear','nearest');
            phiy = F(Ys,Xs,Zs) + Ys;
            F = griddedInterpolant({yI,xI,zI},phiz-ZI,'linear','nearest');
            phiz = F(Ys,Xs,Zs) + Zs;
            
            % find detjac
            [phix_x,phix_y,phix_z] = gradient3d(phix,dxI(1),dxI(2),dxI(3));
            [phiy_x,phiy_y,phiy_z] = gradient3d(phiy,dxI(1),dxI(2),dxI(3));
            [phiz_x,phiz_y,phiz_z] = gradient3d(phiz,dxI(1),dxI(2),dxI(3));
            detjac = phix_x.*(phiy_y.*phiz_z - phiy_z.*phiz_y) ...
                - phix_y.*(phiy_x.*phiz_z - phiy_z.*phiz_x) ...
                + phix_z.*(phiy_x.*phiz_y - phiy_y.*phiz_x);
            
            
            Aphix = A(1,1)*phix + A(1,2)*phiy + A(1,3)*phiz + A(1,4);
            Aphiy = A(2,1)*phix + A(2,2)*phiy + A(2,3)*phiz + A(2,4);
            Aphiz = A(3,1)*phix + A(3,2)*phiy + A(3,3)*phiz + A(3,4);
            
            
            
            % pull back error
            F = griddedInterpolant({yJp,xJp,zJp},errDfp,'linear','nearest');
            lambda = F(Aphiy,Aphix,Aphiz).*detjac.*(-1*det(A)/sigmaM.^2);
            
            % multiply by gradient
            [I_x,I_y,I_z] = gradient3d(It(:,:,:,t),dxI(1),dxI(2),dxI(3));
            grad_x = lambda.*I_x;
            grad_y = lambda.*I_y;
            grad_z = lambda.*I_z;
            
            % smooth
            grad_x = ifftn((vtxhat(:,:,:,t)/sigmaR.^2 + fftn(grad_x).*K).*Kp,'symmetric');
            grad_y = ifftn((vtyhat(:,:,:,t)/sigmaR.^2 + fftn(grad_y).*K).*Kp,'symmetric');
            grad_z = ifftn((vtzhat(:,:,:,t)/sigmaR.^2 + fftn(grad_z).*K).*Kp,'symmetric');
            
            squash = 1;
            if squash
            % update
            stepx = eV*down*grad_x;
            stepy = eV*down*grad_y;
            stepz = eV*down*grad_z;
            
            % squash it            
            normstep = sqrt(stepx.^2 + stepy.^2 + stepz.^2);
            squashnorm = Vmax*normstep./(Vmax + normstep);
            stepx = stepx./normstep.*squashnorm;
            stepy = stepy./normstep.*squashnorm;
            stepz = stepz./normstep.*squashnorm;
            
            
            
            vtx(:,:,:,t) = vtx(:,:,:,t) - stepx;
            vty(:,:,:,t) = vty(:,:,:,t) - stepy;
            vtz(:,:,:,t) = vtz(:,:,:,t) - stepz;
            
            end
            
            if velocity_GN
                % when do I think this is accurate?
                % I think when sigmaR is small (then the correction is
                % small)
                
                % no preconditioner, first thing just smooth and multiply
                % by sigmaR (i.e. apply A^{-1})
                grad_x = sigmaR^2*ifftn((vtxhat(:,:,:,t)/sigmaR.^2 + fftn(grad_x).*K),'symmetric');
                grad_y = sigmaR^2*ifftn((vtyhat(:,:,:,t)/sigmaR.^2 + fftn(grad_y).*K),'symmetric');
                grad_z = sigmaR^2*ifftn((vtzhat(:,:,:,t)/sigmaR.^2 + fftn(grad_z).*K),'symmetric');
                % now multiply by norm of derivative and detjac and weight
                F = griddedInterpolant({yJp,xJp,zJp},WMp,'linear','nearest');
                WMt = F(Aphiy,Aphix,Aphiz);
                factor = (I_x.^2 + I_y.^2 + I_z.^2)/sigmaM^2.*WMt.*detjac;
                grad_x_ = grad_x.*factor;
                grad_y_ = grad_y.*factor;
                grad_z_ = grad_z.*factor;
                % smooth again (don't forget sigmaR)
                grad_x_ = ifftn(fftn(grad_x_).*K,'symmetric')*sigmaR^2;
                grad_y_ = ifftn(fftn(grad_y_).*K,'symmetric')*sigmaR^2;
                grad_z_ = ifftn(fftn(grad_z_).*K,'symmetric')*sigmaR^2;
                % subract
                grad_x = grad_x - grad_x_;
                grad_y = grad_y - grad_y_;
                grad_z = grad_z - grad_z_;
                % update
                vtx(:,:,:,t) = vtx(:,:,:,t) - grad_x*eV;
                vty(:,:,:,t) = vty(:,:,:,t) - grad_y*eV;
                vtz(:,:,:,t) = vtz(:,:,:,t) - grad_z*eV;

            end
                        
        end
        
        
                
        % update
        Ai(1:3,1:4) = Ai(1:3,1:4) - eA * Aistep;
        A = inv(Ai);
        
        % update contrast
        PI = [sum(WM(:)),sum(WA(:)),sum(WB(:))]/numel(WM);
        muB = sum(J(:).*WB(:))/sum(WB(:));
        muA = sum(J(:).*WA(:))/sum(WA(:));
        

        drawnow;
        
    end % end of gradient descent iteration loop it
    
    
    
end % of downsampling loop (downloop)

A = inv(Ai);