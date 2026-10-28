function c = PGA_GeoProd(x,y)
%PGA_GEOPROD Product in the 3D Euclidean PGA even subalgebra.
% Ordered basis: [1,e23,e31,e12,e01,e02,e03,e0123]. e0^2=0.
% Adapted utility; SPDX-License-Identifier: Apache-2.0.
% Modified for the support sample: conventions and input validation clarified;
% row/column inputs normalized. Original scalar product formulas retained.
validateattributes(x,{'numeric'},{'vector','real','finite'});
validateattributes(y,{'numeric'},{'vector','real','finite'});
x=double(x(:)); y=double(y(:));
if numel(x)==6, x=[0;x;0]; end
if numel(y)==6, y=[0;y;0]; end
if numel(x)~=8 || numel(y)~=8
    error('PGA:Size','Inputs must contain six or eight coefficients.');
end
c =[x(1)*y(1) - x(2)*y(2) - x(3)*y(3) - x(4)*y(4)
    x(1)*y(2) + x(2)*y(1) - x(3)*y(4) + x(4)*y(3)
    x(1)*y(3) + x(3)*y(1) + x(2)*y(4) - x(4)*y(2)
    x(1)*y(4) - x(2)*y(3) + x(3)*y(2) + x(4)*y(1)
    x(1)*y(5) + x(5)*y(1) - x(2)*y(8) - x(3)*y(7) + x(4)*y(6) - x(6)*y(4) + x(7)*y(3) - x(8)*y(2)
    x(1)*y(6) + x(6)*y(1) + x(2)*y(7) - x(4)*y(5) + x(5)*y(4) - x(7)*y(2) - x(3)*y(8) - x(8)*y(3)
    x(1)*y(7) - x(2)*y(6) + x(3)*y(5) - x(5)*y(3) + x(6)*y(2) + x(7)*y(1) - x(4)*y(8) - x(8)*y(4)
    x(2)*y(5) + x(5)*y(2) + x(1)*y(8) + x(3)*y(6) + x(6)*y(3) + x(8)*y(1) + x(4)*y(7) + x(7)*y(4)];
end
