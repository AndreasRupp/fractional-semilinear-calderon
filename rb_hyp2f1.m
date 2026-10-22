function f = rb_hyp2f1(a,b,c,z)
if any(z(:)<0 | z(:)>1)
    error('rb_hyp2f1:argument','The argument must lie in [0,1].');
end
f = zeros(size(z));
lo = z<=0.5;
f(lo) = series(a,b,c,z(lo));
hi = ~lo;
if any(hi(:))
    t = 1-z(hi); k = c-a-b;
    if abs(k-round(k))<1e-12
        error('rb_hyp2f1:parameters','An integer continuation exponent is unsupported.');
    end
    A = gamma(c)*gamma(k)/(gamma(c-a)*gamma(c-b));
    B = gamma(c)*gamma(-k)/(gamma(a)*gamma(b));
    f(hi) = A*series(a,b,1-k,t)+B*t.^k.*series(c-a,c-b,1+k,t);
end
end
function y = series(a,b,c,z)
y = ones(size(z)); t = y;
for j=1:200
    t = t.*z*((a+j-1)*(b+j-1)/((c+j-1)*j));
    y = y+t;
    if all(abs(t(:))<=2e-16*max(abs(y(:)),1))
        return
    end
end
error('rb_hyp2f1:series','The hypergeometric series did not converge.');
end
