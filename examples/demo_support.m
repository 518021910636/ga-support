function demo_support
%DEMO_SUPPORT Toy algebra and incidence examples; no application model.

one = [1; zeros(7,1)];
e23 = zeros(8,1); e23(2) = 1;
product = PGA_GeoProd(e23,e23);
assert(norm(product + one) < 1e-12);
fprintf('PGA: e23 * e23 = %.0f\n', product(1));

p = cga5('point',[1,2,3]);
q = cga5('point',[-1,0,2]);
inner = cga5('lc',p,q);
fprintf('CGA: squared Euclidean distance = %.0f\n', -2*inner(1));

sphere = cga5('sphere',[0,0,0],2);
surfacePoint = cga5('point',[2,0,0]);
residual = cga5('lc',surfacePoint,sphere);
fprintf('CGA: sphere-incidence residual = %.0f\n', residual(1));

roundTrip = cga5('down',p);
fprintf('CGA: point round trip = [%g %g %g]\n', roundTrip);
end
