function u0=sj_harmonic(r,reference)
r=r(:);y=reference.inputNodes;w=reference.inputWeights;f=reference.inputValues;s=reference.s;
factor=(2*sin(pi*s)/pi)*((1./bsxfun(@minus,y'.^2,r.^2))*(w.*y.*f./(y.^2-1).^s));
u0=max(0,1-r.^2).^s.*factor;
end
