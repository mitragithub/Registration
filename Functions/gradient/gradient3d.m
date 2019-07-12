function [I_x, I_y, I_z] = gradient3d(I,dx,dy,dz,zero_boundary)
% because a significant amount of time in the gradient function is just
% parsing the input arguments, I will reimplement it
if nargin == 4
    zero_boundary = 0;
end
[ny,nx,nz] = size(I);
I_x = zeros(ny,nx,nz);
I_y = zeros(ny,nx,nz);
I_z = zeros(ny,nx,nz);

% centered difference in middle
I_x(:,2:(nx-1),:) = (I(:,3:nx,:) - I(:,1:(nx-2),:))/(2.0*dx);
I_y(2:(ny-1),:,:) = (I(3:ny,:,:) - I(1:(ny-2),:,:))/(2.0*dy);
I_z(:,:,2:(nz-1)) = (I(:,:,3:nz) - I(:,:,1:(nz-2)))/(2.0*dz);

if ~zero_boundary
% forward difference at start
I_x(:,1,:) = (I(:,2,:) - I(:,1,:))/dx;
I_y(1,:,:) = (I(2,:,:) - I(1,:,:))/dy;
I_z(:,:,1) = (I(:,:,2) - I(:,:,1))/dz;

% backward difference at end
I_x(:,nx,:) = (I(:,nx,:) - I(:,nx-1,:))/dx;
I_y(ny,:,:) = (I(ny,:,:) - I(ny-1,:,:))/dy;
I_z(:,:,nz) = (I(:,:,nz) - I(:,:,nz-1))/dz;
end