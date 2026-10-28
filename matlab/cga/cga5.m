function out = cga5(op,a,b)
% Minimal Cl(4,1): bitmask blades 0:31 in [e1,e2,e3,e+,e-].
% e0=(e--e+)/2, einf=e++e-. No third-party toolbox.
% Adapted support sample: input validation and row normalization added.
% Multivectors always return as rows. No application-specific solver.
persistent idx sg grades
if isempty(idx)
    idx=zeros(32); sg=ones(32); grades=zeros(1,32);
    for i=0:31
        grades(i+1)=sum(bitget(i,1:5));
        for j=0:31
            idx(i+1,j+1)=double(bitxor(i,j))+1; n=0;
            for k=1:5
                if bitget(i,k), n=n+sum(bitget(j,1:k-1)); end
            end
            sg(i+1,j+1)=(-1)^(n+bitget(bitand(i,j),5));
        end
    end
end
switch lower(op)
    case {'gp','wedge','lc'}
        a=validateMV(a); b=validateMV(b);
    case {'reverse','dual','down'}
        a=validateMV(a);
    case 'grade'
        a=validateMV(a);
        validateattributes(b,{'numeric'},{'scalar','integer','>=',0,'<=',5});
end
switch lower(op)
    case 'vector'
        validateattributes(a,{'numeric'},{'vector','numel',5,'real','finite'});
        out=zeros(1,32); out([2 3 5 9 17])=double(a(:).');
    case 'inf'
        out=cga5('vector',[0 0 0 1 1]);
    case 'point'
        validateattributes(a,{'numeric'},{'vector','numel',3,'real','finite'});
        a=double(a(:).'); r2=sum(a.^2);
        out=cga5('vector',[a (r2-1)/2 (r2+1)/2]);
    case 'sphere'
        validateattributes(b,{'numeric'},{'scalar','real','finite','nonnegative'});
        out=cga5('point',a)-0.5*double(b)^2*cga5('inf');
    case 'plane'
        validateattributes(a,{'numeric'},{'vector','numel',3,'real','finite'});
        validateattributes(b,{'numeric'},{'scalar','real','finite'});
        if norm(a)==0, error('CGA:Plane','Plane normal must be nonzero'); end
        out=cga5('vector',[double(a(:).') double(b) double(b)]); % n.x=d
    case 'grade'
        out=a.*(grades==b);
    case 'reverse'
        out=a.*((-1).^(grades.*(grades-1)/2));
    case 'dual'
        I=zeros(1,32); I(32)=-1; out=cga5('gp',a,I);
    case {'gp','wedge','lc'}
        out=zeros(1,32);
        ia=find(a); ib=find(b);
        for i=ia
            for j=ib
                k=idx(i,j);
                if strcmpi(op,'wedge') && bitand(i-1,j-1)~=0, continue; end
                if strcmpi(op,'lc') && (grades(i)>grades(j) || grades(k)~=grades(j)-grades(i)), continue; end
                out(k)=out(k)+a(i)*b(j)*sg(i,j);
            end
        end
    case 'down'
        w=a(17)-a(9);
        if abs(w)<1e-14, error('CGA:Infinity','Cannot dehomogenize infinite point'); end
        out=a([2 3 5])/w;
    otherwise
        error('CGA:Operation','Unknown operation');
end
end

function a=validateMV(a)
validateattributes(a,{'numeric'},{'vector','numel',32,'real','finite'});
a=double(a(:).');
end
