function [g,g0,mu]=sj_flux(rho,reference)
r=reference.integrationNodes;w=reference.integrationWeights;P=reference.integrationPolynomials;s=reference.s;d=reference.d;
kernel=rb_kernel(rho(:),r,d,s);correction=P*reference.coefficients;u0=sj_harmonic(r,reference);weight=(1-r.^2).^s;factor0=u0./weight;
mu=-kernel*(w.*correction);g0=-kernel*(w.*factor0);
y=reference.inputNodes;wy=reference.inputWeights;fy=reference.inputValues;
g0=g0-rb_kernel(rho(:),y,d,s)*(wy.*fy.*y.^(d-1));g=g0+mu;
end
