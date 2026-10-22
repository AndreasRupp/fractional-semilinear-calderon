function report=check_implementation(p)
if nargin<1,p=2;end
assert(isscalar(p)&&ismember(p,[1,2,3]),'Choose p = 1, 2 or 3.');
root=fileparts(mfilename('fullpath'));addpath(root);
outdir=fullfile(root,'verification');if ~exist(outdir,'dir'),mkdir(outdir);end
logfile=fullfile(outdir,'implementation_check.txt');
if p~=2,logfile=fullfile(outdir,sprintf('implementation_check_p%d.txt',p));end
if exist(logfile,'file'),delete(logfile);end
diary(logfile);cleanup=onCleanup(@()diary('off'));
fprintf('Saved-result implementation checks for p = %d; runtime %s\n',p,version);
check_radial_core;
reportCells=cell(2,1);noiseCount=0;runCount=0;
for d=1:2
 folder=fullfile(root,'output',sprintf('example%d',d));
 if p~=2,folder=fullfile(folder,sprintf('p%d',p));end
 saved=load(fullfile(folder,'results.mat'),'R');R=saved.R;cfg=rb_config(d,p);
 assert(isequal(R.cfg,cfg),'Saved configuration differs from rb_config.');
 assert(numel(R.runs)==55,'Expected 30 noise runs and 25 mesh runs.');
 assert(isequal(cfg.mesh_n,[25,50,100,200,400]));
 assert(isequal(cfg.noise_levels,[1e-1,1e-3,1e-5,1e-7,1e-9,1e-10]));
 assert(isequal(cfg.seeds,20261012:20261016));
 assert(cfg.noise_h==.0025 && cfg.mesh_delta==1e-10);
 assert(cfg.s==.35 && cfg.s<.5 && cfg.p==p);
 assert(cfg.q_radius==.75 && cfg.D==.8 && cfg.q_amplitude==1);
 assert(R.reference.degree==256 && R.reference.converged && R.reference.positiveBranch);
 assert(size(R.data.projection,1)==200 && size(R.data.projection,2)==1200);
 sphere=2*pi^(d/2)/gamma(d/2);expectedMass=[];
 for j=1:size(cfg.observation_intervals,1)
  e=linspace(cfg.observation_intervals(j,1),cfg.observation_intervals(j,2),cfg.observation_intervals(j,3)+1)';
  expectedMass=[expectedMass;sphere/d*diff(e.^d)];
 end
 assert(norm(R.data.mass-expectedMass,inf)<2e-14,'Observation masses must use physical radial volume.');
 assert(abs(sqrt(sum(R.data.mass.*R.data.g.^2))-R.data.flux_norm)<1e-13);
 assert(norm(R.data.g-R.data.g0-R.data.mu)<1e-12,'Full/background/correction responses are inconsistent.');
 isnoise=strcmp({R.runs.study},'noise');ismesh=strcmp({R.runs.study},'mesh');
 assert(sum(isnoise)==30 && sum(ismesh)==25);
 noiseRuns=R.runs(isnoise);
 assert(all([noiseRuns.discrepancy_attained]),'A noise-study discrepancy failed.');
 assert(~any([noiseRuns.alpha_at_lower_endpoint]),'A noise run reached the alpha search floor.');
 assert(all([noiseRuns.residual] <= [noiseRuns.target]*(1+1e-6)));
 assert(all(strcmp({noiseRuns.parameter_rule},'discrepancy')));
 maxErrorDifference=0;maxNoiseRelativeDefect=0;maxResidualAbsoluteDefect=0;
 maxPostResidual=0;maxCellFitDefect=0;
 scalarNames={'alpha','tau','beta','visibility','relative_q_error','relative_u_error', ...
  'relative_w_error','residual','target','noise_norm','data_norm'};
 for n=cfg.mesh_n
  idx=find([R.runs.n]==n);z=R.runs(idx);
  [K,M,edges]=rb_energy(n,d,cfg.s);
  B=R.data.projection*rb_observation(R.data.rho,edges,d,cfg.s);
  nd=round(cfg.D*n);ed=edges(1:nd+1);
  [t,w]=rb_gauss(48,0,1);rr=ed(1:end-1)+diff(ed)*t';
  ww=sphere*diff(ed)*w'.*rr.^(d-1);
  qt=zeros(size(rr));inside=rr<cfg.q_radius;
  qt(inside)=cfg.q_amplitude*exp(1-1./(1-(rr(inside)/cfg.q_radius).^2));
  qnorm2=sum(sum(ww.*qt.^2));
  [tf,wf]=rb_gauss(cfg.coefficient_order,0,1);rf=ed(1:end-1)+diff(ed)*tf';
  wf=sphere*diff(ed)*wf'.*rf.^(d-1);vol=sum(wf,2);
  u0=reshape(sj_harmonic(rf(:),R.reference),size(rf));
  for k=1:numel(z)
   a=z(k);
   for j=1:numel(scalarNames)
    v=a.(scalarNames{j});assert(isscalar(v)&&isfinite(v)&&v>=0,'Nonfinite or negative saved diagnostic.');
   end
   assert(all(isfinite([a.q(:);a.w(:);a.correction(:);a.noise(:)])));
   assert(a.h==1/n && isequal(a.edges,ed));
   normNoise=sqrt(sum(R.data.mass.*a.noise.^2));
   defect=abs(normNoise/(a.delta*R.data.flux_norm)-1);
   maxNoiseRelativeDefect=max(maxNoiseRelativeDefect,defect);
   assert(defect<2e-14,'Noise normalization is inconsistent with delta.');
   assert(abs(a.noise_norm-a.delta*a.data_norm)<1e-25);
   res=sqrt(sum(R.data.mass.*(B*a.correction-R.data.mu-a.noise).^2));
   maxResidualAbsoluteDefect=max(maxResidualAbsoluteDefect,abs(res-a.residual));
   assert(abs(res-a.residual)<1e-13,'Saved weighted residual differs from direct evaluation.');
   post=(M+a.tau*K)*a.w-K*a.correction;
   postdef=norm(post)/max(norm(K*a.correction),eps);
   maxPostResidual=max(maxPostResidual,postdef);assert(postdef<1e-11);
   U=min(cfg.A,max(-cfg.A,u0+a.correction(1:nd)));
   mom=sum(wf.*U.^cfg.p,2);den=sum(wf.*U.^(2*cfg.p),2);
   fitted=-a.w(1:nd).*mom./(den+a.beta*vol);
   fitdef=norm(fitted-a.q)/max(norm(a.q),eps);
   maxCellFitDefect=max(maxCellFitDefect,fitdef);assert(fitdef<1e-12);
   assert(abs(min(den./vol)-a.visibility)<1e-13);
   eq=sqrt(sum(sum(ww.*(a.q-qt).^2))/qnorm2);
   maxErrorDifference=max(maxErrorDifference,abs(eq-a.relative_q_error));
   assert(abs(eq-a.relative_q_error)<1e-7,'Coefficient error changes materially on refined quadrature.');
  end
 end
 for seed=cfg.seeds
  oneNoise=R.runs(isnoise & [R.runs.seed]==seed);
  oneMesh=R.runs(ismesh & [R.runs.seed]==seed);
  assert(isequal([oneNoise.delta],cfg.noise_levels));
  assert(all([oneNoise.h]==cfg.noise_h));
  assert(isequal([oneMesh.n],cfg.mesh_n));
  anchor=oneNoise([oneNoise.delta]==cfg.mesh_delta);assert(numel(anchor)==1);
  for k=1:numel(oneNoise)
   assert(norm(oneNoise(k).noise/oneNoise(k).delta-anchor.noise/anchor.delta)<1e-12);
   expectedBeta=cfg.beta_scale*mean_visibility_integral(oneNoise(k),R.reference,cfg);
   assert(abs(oneNoise(k).beta-expectedBeta)<1e-22,'Beta does not follow the state-based rule.');
  end
  for k=1:numel(oneMesh)
   a=oneMesh(k);
   assert(a.delta==cfg.mesh_delta && a.alpha==anchor.alpha && a.tau==anchor.tau && a.beta==anchor.beta);
   assert(isequal(a.noise,anchor.noise));
   assert(strcmp(a.parameter_rule,'fixed_finest_mesh'));
  end
  assert(norm(oneMesh(end).q-anchor.q)==0,'Finest mesh run differs from its noise-study anchor.');
 end
 assert(all(abs([noiseRuns.tau]/(cfg.tau_scale*cfg.noise_h)-1)<1e-14));
 csv_check(fullfile(folder,'runs.csv'),R);
 row=struct('dimension',d,'power',p,'runs',numel(R.runs),'noiseRuns',sum(isnoise), ...
  'maxNoiseRelativeDefect',maxNoiseRelativeDefect,'maxResidualAbsoluteDefect',maxResidualAbsoluteDefect, ...
  'maxPostResidual',maxPostResidual,'maxCellFitDefect',maxCellFitDefect, ...
  'maxCoefficientErrorQuadratureDifference',maxErrorDifference, ...
  'minimumSelectedAlpha',min([noiseRuns.alpha]),'maximumSelectedAlpha',max([noiseRuns.alpha]));
 reportCells{d}=row;
 fprintf('\nDimension %d, p = %d: 55 runs, 30 attained noise discrepancies, no alpha-floor hits.\n',d,p);
 fprintf('Exact noise and mesh grids, five paired seeds, physical volume weights: pass.\n');
 fprintf('Finest-grid alpha/tau/beta fixed across mesh widths: pass.\n');
 fprintf('All CSV entries agree with the saved arrays: pass.\n');
 fprintf('Noise-norm relative defect: %.3g; direct-residual defect: %.3g.\n',maxNoiseRelativeDefect,maxResidualAbsoluteDefect);
 fprintf('Weak-Laplacian equation residual: %.3g; coefficient-fit defect: %.3g.\n',maxPostResidual,maxCellFitDefect);
 fprintf('12-to-48-point coefficient-error change: %.3g.\n',maxErrorDifference);
 fprintf('Selected noise-study alpha range: %.8g to %.8g.\n',row.minimumSelectedAlpha,row.maximumSelectedAlpha);
 fprintf('delta, mean relative coefficient error\n');
 for delta=cfg.noise_levels
  select=R.runs(isnoise & [R.runs.delta]==delta);
  fprintf('%.1e, %.10g\n',delta,mean([select.relative_q_error]));
 end
 fprintf('h, mean relative coefficient error (fixed regularization)\n');
 for n=cfg.mesh_n
  select=R.runs(ismesh & [R.runs.n]==n);
  fprintf('%.4g, %.10g\n',1/n,mean([select.relative_q_error]));
 end
 noiseCount=noiseCount+sum(isnoise);runCount=runCount+numel(R.runs);
end
report=vertcat(reportCells{:});
fprintf('\nPASS for p = %d: %d saved runs; all %d noise-study discrepancies attained.\n',p,runCount,noiseCount);
fprintf('Source audit: rb_inverse selects alpha from the measured residual, tau from h, and beta from the reconstructed state.\n');
fprintf('True q, u and w are used only to generate synthetic data and report diagnostic errors.\n');
fprintf('These checks do not prove the coupled asymptotic parameter conditions or guarantee general nonsymmetric recovery.\n');
end

function value=mean_visibility_integral(a,raw,cfg)
n=a.n;nd=numel(a.q);e=a.edges;[t,w]=rb_gauss(cfg.coefficient_order,0,1);
r=e(1:end-1)+diff(e)*t';sphere=2*pi^(cfg.d/2)/gamma(cfg.d/2);
weights=sphere*diff(e)*w'.*r.^(cfg.d-1);
U=reshape(sj_harmonic(r(:),raw),size(r))+a.correction(1:nd);
U=min(cfg.A,max(-cfg.A,U));
value=sum(sum(weights.*U.^(2*cfg.p)))/sum(weights(:));
assert(n==round(1/a.h));
end

function csv_check(path,R)
fid=fopen(path,'r');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
header=fgetl(fid);assert(ischar(header)&&~isempty(strfind(header,'dimension,study,seed')));
for k=1:numel(R.runs)
 line=fgetl(fid);assert(ischar(line));tokens=strsplit(line,',');a=R.runs(k);
 assert(numel(tokens)==16 && strcmp(tokens{2},a.study));
 columns=[1,3:16];actual=zeros(size(columns));
 for j=1:numel(columns),actual(j)=str2double(tokens{columns(j)});end
 expected=[R.d,a.seed,a.h,a.delta,a.alpha,a.tau,a.beta,a.relative_q_error,a.relative_u_error,a.relative_w_error,a.visibility,a.residual,a.noise_norm,a.target,a.discrepancy_attained];
 assert(all(abs(actual-expected)<=5e-15*max(abs(expected),realmin)),'CSV value differs from saved result.');
end
assert(~ischar(fgetl(fid)),'Unexpected additional CSV rows.');
end
