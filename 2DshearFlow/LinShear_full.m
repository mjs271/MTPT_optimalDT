function nu = LinShear_full(x1,x2,D1,D2,al,t)
%LinShear_full Full anisotropic solution from Bolster et al. hypermixing

k=zeros(2,2);
k(1,1)=2*D1.*t + (2/3)*D2*(al^2).*(t.^3);
k(2,1)=D2*al.*(t.^2);
k(1,2)=D2*al.*(t.^2);
k(2,2)=2*D2.*t;

x=[x1; x2];

nu=(exp(-0.5*transpose(x)*(k^-1)*x))/(2*pi.*sqrt(det(k)));
end

