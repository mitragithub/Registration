function W = build_boundary_weight(n,b)
% n is size
% b is size of boundary
if nargin == 0
    n = [59,45];
    b = 10;
end

Wx = ones(n);
Wx(:,1:b) = repmat( 0.5*(1-cos((0:b-1)*pi/b))  , [n(1),1]);
Wx = Wx .* flip(Wx,2);

Wy = ones(n);
Wy(1:b,:) = repmat( 0.5*(1-cos((0:b-1)*pi/b)')  , [1,n(2)]);
Wy = Wy .* flip(Wy,1);


W = Wx.*Wy;