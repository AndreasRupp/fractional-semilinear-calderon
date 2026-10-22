function reference=sj_reference(d,s,n,options)
if nargin<4,options=struct;end
A=setting(options,'amplitude',32);radius=setting(options,'radius',.65);p=setting(options,'power',2);qAmplitude=setting(options,'qAmplitude',1);
inputCentre=setting(options,'inputCentre',1.6);inputRadius=setting(options,'inputRadius',.2);qorder=setting(options,'quadratureOrder',max(400,4*n));
[t,wt]=sj_quadrature(setting(options,'inputOrder',128),0,0);inputNodes=inputCentre+inputRadius*t;inputWeights=inputRadius*wt;inputValues=A*sj_bump(t);
reference=struct('d',d,'s',s,'degree',n,'amplitude',A,'radius',radius,'p',p,'qAmplitude',qAmplitude,'inputCentre',inputCentre,'inputRadius',inputRadius, ...
 'inputNodes',inputNodes,'inputWeights',inputWeights,'inputValues',inputValues,'quadratureOrder',qorder);
[t,wt]=sj_quadrature(qorder,0,0);r=radius*(t+1)/2;sphere=2*pi^(d/2)/gamma(d/2);w=sphere*radius/2*wt.*r.^(d-1);
[P,lambda]=sj_basis(r,n,d,s);Phi=bsxfun(@times,(1-r.^2).^s,P);u0=sj_harmonic(r,reference);q=qAmplitude*sj_bump(r/radius);c=zeros(n+1,1);
residuals=zeros(40,1);converged=false;newtonUpdates=0;stepFactors=zeros(40,1);backtrackEvaluations=zeros(40,1);
for k=1:40
 u=u0+Phi*c;res=lambda.*c+Phi'*(w.*q.*max(u,0).^p);residuals(k)=norm(res./sqrt(lambda));
 if residuals(k)<1e-13,converged=true;break;end
 J=diag(lambda)+Phi'*bsxfun(@times,w.*q.*p.*max(u,0).^(p-1),Phi);step=-(J\res);factor=1;
 for bt=1:30
  trialFactor=factor;trial=c+factor*step;ut=u0+Phi*trial;rr=lambda.*trial+Phi'*(w.*q.*max(ut,0).^p);
  if norm(rr./sqrt(lambda))<=(1-1e-4*factor)*residuals(k),break;end;factor=factor/2;
 end
 c=trial;newtonUpdates=newtonUpdates+1;stepFactors(newtonUpdates)=trialFactor;backtrackEvaluations(newtonUpdates)=bt;
end
reference.coefficients=c;reference.lambda=lambda;reference.residual=residuals(k);reference.residualHistory=residuals(1:k);reference.converged=converged;
reference.newtonUpdates=newtonUpdates;reference.jacobianSolves=newtonUpdates;reference.stepFactors=stepFactors(1:newtonUpdates);reference.backtrackEvaluations=backtrackEvaluations(1:newtonUpdates);
[rx,wx]=sj_quadrature(qorder,s,d/2-1);rx=sqrt((rx+1)/2);wx=2^(-s-d/2-1)*wx;
reference.integrationNodes=rx;reference.integrationWeights=wx;reference.integrationPolynomials=sj_basis(rx,n,d,s);
rr=linspace(0,1,2001)';u=sj_harmonic(rr,reference)+(1-rr.^2).^s.*(sj_basis(rr,n,d,s)*c);
reference.minimumState=min(u);reference.maximumState=max(u);reference.positiveBranch=reference.minimumState>=-1e-10&&reference.maximumState<=A+1e-10;
assert(converged&&reference.positiveBranch,'Radial reference failed forward residual or state-bound checks.');
end
function value=setting(o,name,default)
if isfield(o,name),value=o.(name);else,value=default;end
end
