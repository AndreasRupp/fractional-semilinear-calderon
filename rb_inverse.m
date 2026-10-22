function out = rb_inverse(model,delta,seed)
if nargin<3,seed=20261012;end
rng(seed,'twister');data=model.data;xi=randn(size(data.g));
noise=delta*data.flux_norm*xi/sqrt(sum(data.mass.*xi.^2));
y=sqrt(data.mass).*(data.mu+noise);z=model.Us'*y;
yperp=y-model.Us*z;s=model.sigma;
res2=sum((bsxfun(@times,z,bsxfun(@rdivide,model.alphas,s.^2+model.alphas))).^2,1)+sum(yperp.^2);
res=sqrt(res2);target=model.discrepancy_factor*delta*data.flux_norm;
i=find(res<=target,1,'first');attained=~isempty(i);
if ~attained,i=numel(model.alphas);end
alpha=model.alphas(i);r=model.R\(model.Vs*((s./(s.^2+alpha)).*z));
w=model.post\(model.post'\(model.K*r));
U=model.u0+r(1:model.nd);U=min(model.amplitude,max(-model.amplitude,U));
mom=sum(model.weights.*U.^model.p,2);den=sum(model.weights.*U.^(2*model.p),2);
vol=sum(model.weights,2);visibility=min(den./vol);beta=model.beta_scale*sum(den)/sum(vol);
if ~isempty(model.fixed_beta),beta=model.fixed_beta;end
q=-w(1:model.nd).*mom./(den+beta*vol);
qerror=sqrt(sum(sum(model.weights.*(q-model.q).^2))/sum(sum(model.weights.*model.q.^2)));
stateerror=sqrt(sum(sum(model.weights.*(U-model.u).^2))/sum(sum(model.weights.*model.u.^2)));
posterror=sqrt(sum(sum(model.weights.*(w(1:model.nd)-model.wtrue).^2))/sum(sum(model.weights.*model.wtrue.^2)));
actual=sqrt(sum(data.mass.*(model.B*r-data.mu-noise).^2));
if numel(model.alphas)==1,rule='fixed';else,rule='Morozov';end
out=struct('d',model.d,'n',model.n,'h',model.h,'delta',delta,'seed',seed,'alpha',alpha,'tau',model.tau,'beta',beta,'visibility',visibility, ...
 'parameter_rule',rule,'alpha_at_lower_endpoint',numel(model.alphas)>1&&i==numel(model.alphas),'relative_q_error',qerror,'relative_u_error',stateerror,'relative_w_error',posterror,'residual',actual,'target',target,'discrepancy_attained',attained, ...
 'q',q,'w',w,'correction',r,'edges',model.edges(1:model.nd+1),'noise',noise,'noise_norm',delta*data.flux_norm,'data_norm',data.flux_norm);
end
