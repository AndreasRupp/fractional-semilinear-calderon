function [u,u0,q,w]=sj_evaluate(r,reference)
r=r(:);u0=sj_harmonic(r,reference);P=sj_basis(r,reference.degree,reference.d,reference.s);u=u0+max(0,1-r.^2).^reference.s.*(P*reference.coefficients);
q=reference.qAmplitude*sj_bump(r/reference.radius);w=-q.*u.^reference.p;
end
