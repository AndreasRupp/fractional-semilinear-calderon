function [x,w]=sj_quadrature(n,a,b)
j=(0:n-1)';den=(2*j+a+b).*(2*j+a+b+2);d=(b*b-a*a)./den;if abs(a+b)<1e-14,d(1)=(b-a)/(a+b+2);end
k=(1:n-1)';e=2./(2*k+a+b).*sqrt(k.*(k+a).*(k+b).*(k+a+b)./((2*k+a+b-1).*(2*k+a+b+1)));
[V,D]=eig(diag(d)+diag(e,1)+diag(e,-1));[x,p]=sort(diag(D));w=exp((a+b+1)*log(2)+gammaln(a+1)+gammaln(b+1)-gammaln(a+b+2))*V(1,p)'.^2;
end
