% gaussian mixture model

function [L,mu] = kmeans(I,OPT)
if nargin == 1
    OPT = struct;
end
%default options
image = 0; % show an image, dimensions are N1xN2xC channels, versus NxC
niter = 10;
M = 2; % components
draw = 37;

% update options
if isfield(OPT,'image')
    image = OPT.image;
end
if isfield(OPT,'niter')
    niter = OPT.niter;
end
if isfield(OPT,'M')
    M = OPT.M;
end
if isfield(OPT,'draw')
    draw = OPT.draw;
end


if image
    imagesize = size(I);
    N = imagesize(1)*imagesize(2);
    C = size(I,3);
    I = reshape(I,N,C);
else
    N = size(I,1);
    C = size(I,2);
end






% initialize parameters
% I will rearrange my data as follows
% first axis realizations
% second axis vector components
% third axis mixture component
% or for covariance, the first and second are vector components

mu = zeros(1,C,M);
for c = 1 : C
    qs = linspace(min(I(:,c)),max(I(:,c)),M+2);
    mu(1,c,:) = qs(2:end-1);
end
if isfield(OPT,'mu')
    mu = OPT.mu;
end







if image
    if draw
        h = danfigure(draw);
%         h = danfigure(37);
        subplot(1,2,1)
        imagesc(reshape(I,imagesize))
        axis image
    end
end



for it = 1 : niter
    % distance
    I0 = bsxfun(@minus,I,mu);
    d2 = sum(I0.^2,2);
    
    % assign
    [mind2,minind] = min(d2,[],3);
    oh = zeros(N,1,M);
    for m = 1 : M
        oh(:,:,m) = minind==m;
    end
    
    
    
    
    
    % update params, weight by psample
    mu = bsxfun(@rdivide, sum(bsxfun(@times,I,oh),1)  , sum(oh,1));

    
    
    
    if image
        if draw
            danfigure(h);
            subplot(1,2,2);
            if M == 2
                imagesc(reshape(cat(3,oh(:,:,1),oh(:,:,2),oh(:,:,1)),imagesize(1),imagesize(2),3))
            else
                imagesc(reshape(oh(:,:,1:3),imagesize(1),imagesize(2),3))
            end
            axis image
            drawnow
        end
    end
    
    
end

L = oh;
if image
    L = reshape(oh,imagesize(1),imagesize(2),M);
end