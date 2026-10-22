function check_radial_core
s=.35; C=4^s*gamma(.5+s)/(sqrt(pi)*abs(gamma(-s)));
[R,T]=ndgrid([.02,.2,.7,1],[.03,.3,.7,1]);
E=rb_ball_energy(R(:,1),T(1,:),1,s);
beta=1-2*s;
exact=C/(s*beta)*((R+T).^beta-abs(R-T).^beta);
assert(norm(E-exact,'fro')/norm(exact,'fro')<2e-13);
for d=1:2
    for n=[25,100,400]
        [K,M,e]=rb_energy(n,d,s); chol(K);
        assert(max(max(K-diag(diag(K))))<2e-12);
        assert(min(sum(K,2))>0);
        rho=[1.01;1.1;1.3;2]; B=rb_observation(rho,e,d,s);
        assert(max(B(:))<0);
        if d==1
            F=@(r) -C/(2*s)*((rho-r).^(-2*s)-(rho+r).^(-2*s));
            Bexact=diff(F(e'),1,2);
            assert(norm(B-Bexact,'fro')/norm(Bexact,'fro')<2e-12);
        end
        fprintf('d=%d n=%d: positive stiffness, sign checks and exact kernel checks pass\n',d,n);
    end
end
[r,w]=rb_gauss(100,0,.7); rho=[1.01;1.5;2];
for d=1:2
    I=-rb_kernel(rho,r,d,s)*(w.*r.^(d-1));
    B=rb_observation(rho,[0;.7],d,s);
    assert(norm(I-B)/norm(B)<2e-13);
end
fprintf('Angular kernel integration agrees with exact ball observations.\n');
end
