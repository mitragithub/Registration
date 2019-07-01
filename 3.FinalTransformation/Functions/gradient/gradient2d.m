function [I_x, I_y] = gradient2d(I,dx,dy)
% because a significant amount of time in the gradient function is just
% parsing the input arguments, I will reimplement it
[ny,nx] = size(I);
I_x = zeros(ny,nx);
I_y = zeros(ny,nx);

% centered difference in middle
I_x(:,2:(nx-1)) = (I(:,3:nx) - I(:,1:(nx-2)))/(2.0*dx);
I_y(2:(ny-1),:) = (I(3:ny,:) - I(1:(ny-2),:))/(2.0*dy);

% forward difference at start
I_x(:,1) = (I(:,2) - I(:,1))/dx;
I_y(1,:) = (I(2,:) - I(1,:))/dy;

% backward difference at end
I_x(:,nx) = (I(:,nx) - I(:,nx-1))/dx;
I_y(ny,:) = (I(ny,:) - I(ny-1,:))/dy;