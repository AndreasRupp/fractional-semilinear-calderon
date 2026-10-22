function E = rb_ball_energy(R,T,d,s)
if s<=0 || s>=0.5 || ~ismember(d,[1,2])
    error('rb_ball_energy:parameters','Use d=1 or 2 and 0<s<1/2.');
end
R=R(:); T=T(:)';
a=min(R,T); b=max(R,T); nu=d/2;
E=zeros(size(a)); ind=a>0;
z=(a(ind)./b(ind)).^2;
S=2*pi^nu/gamma(nu);
factor=S*gamma(nu+s)/(2^(1-2*s)*gamma(1-s)*gamma(nu+1));
E(ind)=factor*a(ind).^d.*b(ind).^(-2*s).*rb_hyp2f1(nu+s,s,nu+1,z);
end
