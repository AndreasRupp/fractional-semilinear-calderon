function model = rb_prepare(n,reference,cfg,data)
if nargin<3,cfg=struct;end
if nargin<4,data=rb_data(reference,cfg);end
d=reference.raw.d;s=reference.raw.s;p=reference.raw.p;
[K,M,edges]=rb_energy(n,d,s);mass=diag(M);
B=data.projection*rb_observation(data.rho,edges,d,s);
R=chol(K+M);A=bsxfun(@times,sqrt(data.mass),B)/R;
[Us,Ss,Vs]=svd(A,'econ');rankdim=min(size(A));
Us=Us(:,1:rankdim);Vs=Vs(:,1:rankdim);sigma=diag(Ss(1:rankdim,1:rankdim));
D=getoption(cfg,'D',.8);nd=round(n*D);
assert(abs(nd/n-D)<1e-12,'The reconstruction radius must be fitted to the radial mesh.');
[t,w]=rb_gauss(getoption(cfg,'coefficient_order',12),0,1);
r=edges(1:nd)+(edges(2:nd+1)-edges(1:nd))*t';
weights=2*pi^(d/2)/gamma(d/2)*(edges(2:nd+1)-edges(1:nd))*w'.*r.^(d-1);
[u,u0,q,wt]=sj_evaluate(r(:),reference.raw);
u=reshape(u,size(r));u0=reshape(u0,size(r));q=reshape(q,size(r));wt=reshape(wt,size(r));
tau=getoption(cfg,'tau',getoption(cfg,'tau_scale',1e-4)/n);
post=chol(M+tau*K);
model=struct('n',n,'h',1/n,'d',d,'s',s,'p',p,'D',D,'nd',nd,'K',K,'M',M,'mass',mass,'edges',edges,'B',B,'R',R,'Us',Us,'sigma',sigma,'Vs',Vs, ...
 'r',r,'weights',weights,'u',u,'u0',u0,'q',q,'wtrue',wt,'tau',tau,'post',post,'amplitude',reference.raw.amplitude,'data',data, ...
 'fixed_beta',getoption(cfg,'beta',[]),'beta_scale',getoption(cfg,'beta_scale',1e-10),'alphas',getoption(cfg,'alphas',10.^linspace(0,-16,161)),'discrepancy_factor',getoption(cfg,'discrepancy_factor',1.05));
end
function v=getoption(c,n,default)
if isfield(c,n),v=c.(n);else,v=default;end
end
