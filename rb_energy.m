function [K,M,edges] = rb_energy(n,d,s)
edges=linspace(0,1,n+1)';
E=rb_ball_energy(edges,edges,d,s);
K=diff(diff(E,1,1),1,2);
K=(K+K')/2;
volume=pi^(d/2)/gamma(1+d/2);
M=diag(volume*diff(edges.^d));
end
