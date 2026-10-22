function G = rb_kernel(r,t,d,s)
r=r(:); t=t(:)';
a=min(r,t); b=max(r,t);
if any(a(:)==b(:))
    error('rb_kernel:separation','The radial arguments must be separated.');
end
nu=d/2;
C=4^s*gamma(nu+s)/(pi^nu*abs(gamma(-s)));
S=2*pi^nu/gamma(nu);
G=C*S*b.^(-d-2*s).*rb_hyp2f1(nu+s,1+s,nu,(a./b).^2);
end
