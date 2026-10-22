function B = rb_observation(rho,edges,d,s)
rho=rho(:); edges=edges(:)';
if any(rho<=max(edges))
    error('rb_observation:separation','Observation radii must exceed all cell radii.');
end
nu=d/2;
C=4^s*gamma(nu+s)/(pi^nu*abs(gamma(-s)));
V=pi^nu/gamma(nu+1);
z=(edges./rho).^2;
ball=-C*V*(edges.^d).*rho.^(-d-2*s).*rb_hyp2f1(nu+s,1+s,nu+1,z);
B=diff(ball,1,2);
end
