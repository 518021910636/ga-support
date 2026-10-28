function report = test_support
%TEST_SUPPORT Algebra conventions and geometric identities on synthetic data.
previousRng = rng;
cleanup = onCleanup(@() rng(previousRng)); %#ok<NASGU>
rng(17,'twister');
tol = 1e-10;
one = [1;zeros(7,1)];

% Independent bit-mask reference for all 64 even-basis PGA products.
masks = [0,12,10,6,3,5,9,15];
orientation = [1,1,-1,1,1,1,1,1]; % e31 = -e13
for i = 1:8
    for j = 1:8
        x=zeros(8,1); x(i)=1;
        y=zeros(8,1); y(j)=1;
        expected=zeros(8,1);
        if bitand(bitand(masks(i),masks(j)),1)==0
            swaps=0;
            for bit=0:3
                if bitget(masks(i),bit+1)
                    swaps=swaps+sum(bitget(masks(j),1:bit));
                end
            end
            k=find(masks==bitxor(masks(i),masks(j)));
            expected(k)=(-1)^swaps*orientation(i)*orientation(j)/orientation(k);
        end
        assert(norm(PGA_GeoProd(x,y)-expected)<tol);
    end
end

for k=1:20
    x=randn(8,1); y=randn(8,1); z=randn(8,1);
    assert(norm(PGA_GeoProd(x,one)-x)<tol);
    assert(norm(PGA_Rev(PGA_Rev(x))-x)<tol);
    assert(norm(PGA_GeoProd(PGA_GeoProd(x,y),z)- ...
        PGA_GeoProd(x,PGA_GeoProd(y,z)))<tol);
    assert(norm(PGA_Rev(PGA_GeoProd(x,y))- ...
        PGA_GeoProd(PGA_Rev(y),PGA_Rev(x)))<tol);
    v=randn(1,6);
    assert(norm(PGA_GeoProd(v,v)-PGA_GeoProd([0;v(:);0],[0;v(:);0]))<tol);
    assert(norm(PGA_GeoProd(x.',y.')-PGA_GeoProd(x,y))<tol);
end

% Orthogonal CGA metric and null-basis convention.
metric=[1,1,1,1,-1];
for i=1:5
    a=zeros(1,5); a(i)=1; a=cga5('vector',a);
    for j=1:5
        b=zeros(1,5); b(j)=1; b=cga5('vector',b);
        anti=cga5('gp',a,b)+cga5('gp',b,a);
        expected=zeros(1,32);
        if i==j, expected(1)=2*metric(i); end
        assert(norm(anti-expected)<tol);
    end
end
eo=cga5('point',[0,0,0]); einf=cga5('inf');
dotNull=cga5('lc',eo,einf);
assert(abs(dotNull(1)+1)<tol);
assert(norm(cga5('gp',einf,einf))<tol);

for k=1:20
    x=randn(1,3); y=randn(1,3); center=randn(1,3); radius=0.5+rand;
    p=cga5('point',x); q=cga5('point',y);
    assert(norm(cga5('down',p)-x)<tol);
    assert(norm(cga5('down',-3*p)-x)<tol);
    assert(norm(cga5('gp',p,p))<tol);
    distance=cga5('lc',p,q);
    assert(abs(-2*distance(1)-sum((x-y).^2))<tol);
    sphere=cga5('sphere',center,radius);
    residual=cga5('lc',p,sphere);
    assert(abs(residual(1)-0.5*(radius^2-sum((x-center).^2)))<tol);
    normal=randn(1,3); offset=randn;
    plane=cga5('plane',normal,offset);
    residual=cga5('lc',p,plane);
    assert(abs(residual(1)-(dot(normal,x)-offset))<tol);
    assert(norm(cga5('wedge',p,p))<tol);
    assert(norm(cga5('gp',p.',q.')-cga5('gp',p,q))<tol);
end

for k=1:6
    a=randn(1,32); b=randn(1,32); c=randn(1,32);
    assert(norm(cga5('gp',cga5('gp',a,b),c)- ...
        cga5('gp',a,cga5('gp',b,c)))<tol);
    assert(norm(cga5('reverse',cga5('gp',a,b))- ...
        cga5('gp',cga5('reverse',b),cga5('reverse',a)))<tol);
    assert(norm(cga5('dual',cga5('dual',a))+a)<tol);
    g2=cga5('grade',a,2); g3=cga5('grade',b,3);
    contraction=cga5('lc',g2,g3);
    assert(norm(contraction-cga5('grade',contraction,1))<tol);
    assert(norm(cga5('lc',g3,g2))<tol);
end
assertThrows(@() PGA_GeoProd(zeros(5,1),one));
assertThrows(@() PGA_Rev(zeros(6,1)));
assertThrows(@() cga5('point',[NaN,0,0]));
assertThrows(@() cga5('sphere',[0,0,0],-1));
assertThrows(@() cga5('down',einf));
assertThrows(@() cga5('gp',zeros(1,31),zeros(1,32)));
assertThrows(@() cga5('grade',zeros(1,32),6));

report=struct('passed',true,'pga_basis_products',64, ...
    'pga_random_trials',20,'cga_metric_pairs',25, ...
    'cga_geometry_trials',20,'cga_dense_trials',6,'invalid_input_cases',7);
fprintf('PASS: PGA / CGA support tests\n');
end

function assertThrows(f)
threw=false;
try
    f();
catch
    threw=true;
end
assert(threw,'Expected invalid input to raise an error.');
end
