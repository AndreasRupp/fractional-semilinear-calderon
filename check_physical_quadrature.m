function report=check_physical_quadrature
s=.35;rows=[];
for d=1:2
 pairs=[.3,.7;.99,1.05;.7,2];
 for k=1:size(pairs,1)
  r=pairs(k,1);t=pairs(k,2);a=rb_kernel(r,t,d,s);b=constant(d,s)*physical_angle(r,t,d,d/2+s);
  defect=abs(a-b)/abs(b);assert(defect<1e-9);rows=[rows;d,1,r,t,defect];
 end
 pairs=[1.05,1;1.3,.8;2,.4];
 for k=1:size(pairs,1)
  rho=pairs(k,1);R=pairs(k,2);a=rb_observation(rho,[0;R],d,s);
  b=-constant(d,s)*integral(@(t)arrayfun(@(v)v^(d-1)*physical_angle(rho,v,d,d/2+s),t),0,R,'RelTol',1e-11,'AbsTol',1e-12);
  defect=abs(a-b)/abs(b);assert(defect<1e-9);rows=[rows;d,2,rho,R,defect];
 end
 [t,wt]=sj_quadrature(128,0,0);
 raw=struct('d',d,'s',s,'inputNodes',1.6+.2*t,'inputWeights',.2*wt,'inputValues',16*sj_bump(t));
 for r=[0,.5,.99]
  a=sj_harmonic(r,raw);factor=gamma(d/2)*sin(pi*s)/pi^(d/2+1)*(1-r^2)^s;
  b=factor*integral(@(t)arrayfun(@(v)v^(d-1)*16*sj_bump((v-1.6)/.2)/(v^2-1)^s*physical_angle(r,v,d,d/2),t),1.4,1.8,'RelTol',1e-11,'AbsTol',1e-12);
  defect=abs(a-b)/abs(b);assert(defect<1e-9);rows=[rows;d,3,r,0,defect];
 end
end
pairs=[.2,.7;.7,1;.7,.7;1,1];
for k=1:size(pairs,1)
 a=pairs(k,1);b=pairs(k,2);value=rb_ball_energy(a,b,2,s);
 if a==b
  exponent=1/(1-2*s);
  inside=exponent*(2*a)^(1-2*s)*integral(@(x)deficit_over_distance(2*a*x.^exponent,a),0,1,'RelTol',1e-11,'AbsTol',1e-12);
 else
  inside=integral(@(z)disk_deficit(z,a,b).*z.^(-1-2*s),b-a,a+b,'RelTol',1e-11,'AbsTol',1e-12);
 end
 physical=constant(2,s)*2*pi*(inside+pi*a^2*(a+b)^(-2*s)/(2*s));
 defect=abs(value-physical)/abs(physical);assert(defect<1e-9);rows=[rows;2,4,a,b,defect];
end
report=array2table(rows,'VariableNames',{'dimension','test','radius1','radius2','relative_difference'});
disp(report);
fprintf('Tests 1/2/3/4: physical angular kernel, ball observation, Poisson lift, full-space disk energy.\n');
fprintf('These independent integration checks do not bound the continuum forward error.\n');
end
function value=constant(d,s)
value=4^s*gamma(d/2+s)/(pi^(d/2)*abs(gamma(-s)));
end
function value=physical_angle(r,t,d,power)
if d==1
 value=abs(r-t)^(-2*power)+(r+t)^(-2*power);
else
 value=2*integral(@(theta)((r-t)^2+4*r*t*sin(theta/2).^2).^(-power),0,pi,'RelTol',1e-12,'AbsTol',1e-12);
end
end
function value=deficit_over_distance(z,a)
value=2*a*ones(size(z));idx=z>0;
value(idx)=2*a^2*asin(z(idx)/(2*a))./z(idx)+.5*sqrt(max(0,4*a^2-z(idx).^2));
end
function value=disk_deficit(z,a,b)
ca=min(1,max(-1,(z.^2+a^2-b^2)./(2*z*a)));
cb=min(1,max(-1,(z.^2+b^2-a^2)./(2*z*b)));
area=a^2*acos(ca)+b^2*acos(cb)-.5*sqrt(max(0,(-z+a+b).*(z+a-b).*(z-a+b).*(z+a+b)));
value=pi*a^2-area;
end
