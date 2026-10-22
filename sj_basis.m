function [P,lambda]=sj_basis(r,n,d,s)
a=s;b=d/2-1;z=2*r(:).^2-1;P=zeros(numel(z),n+1);P(:,1)=1;if n>0,P(:,2)=.5*((a-b)+(a+b+2)*z);end
for k=2:n,t=2*k+a+b;P(:,k+1)=((t-1)*((t*(t-2))*z+a*a-b*b).*P(:,k)-2*(k+a-1)*(k+b-1)*t*P(:,k-1))/(2*k*(k+a+b)*(t-2));end
j=(0:n)';sphere=2*pi^(d/2)/gamma(d/2);
logh=(a+b+1)*log(2)-log(2*j+a+b+1)+gammaln(j+a+1)+gammaln(j+b+1)-gammaln(j+1)-gammaln(j+a+b+1);
h=sphere*2^(-s-d/2-1)*exp(logh);P=bsxfun(@rdivide,P,sqrt(h)');
lambda=exp(2*s*log(2)+gammaln(j+s+1)+gammaln(j+s+d/2)-gammaln(j+1)-gammaln(j+d/2));
end
