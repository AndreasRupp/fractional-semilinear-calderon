function data = rb_data(reference,cfg)
if nargin<2,cfg=struct;end
intervals=getoption(cfg,'observation_intervals',[1.05,1.3,120;2,3,80]);
order=getoption(cfg,'observation_order',6);d=reference.raw.d;
a=[];b=[];
for j=1:size(intervals,1)
    e=linspace(intervals(j,1),intervals(j,2),intervals(j,3)+1)';
    a=[a;e(1:end-1)];b=[b;e(2:end)];
end
[t,v]=rb_gauss(order,0,1);r=a+(b-a)*t';
weights=2*pi^(d/2)/gamma(d/2)*(b-a)*v'.*r.^(d-1);
vol=sum(weights,2);projection=weights./vol;
rho=reshape(r',[],1);w=reshape(projection',[],1);
index=repelem((1:numel(a))',order);Q=sparse(index,(1:numel(rho))',w,numel(a),numel(rho));
[g,g0,mu]=sj_flux(rho,reference.raw);
data=struct('rho',rho,'projection',Q,'mass',vol,'g',Q*g,'g0',Q*g0,'mu',Q*mu,'intervals',intervals,'order',order);
data.flux_norm=sqrt(sum(vol.*data.g.^2));
end
function v=getoption(c,n,default)
if isfield(c,n),v=c.(n);else,v=default;end
end
