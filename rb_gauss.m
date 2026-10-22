function [x,w] = rb_gauss(n,a,b)
k=(1:n-1)'; t=k./sqrt(4*k.^2-1);
[V,D]=eig(diag(t,1)+diag(t,-1));
[x,ind]=sort(diag(D)); w=2*V(1,ind)'.^2;
x=a+(b-a)*(x+1)/2; w=(b-a)*w/2;
end
