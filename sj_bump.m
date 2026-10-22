function q=sj_bump(r)
q=zeros(size(r));inside=abs(r)<1;q(inside)=exp(1-1./(1-r(inside).^2));
end
