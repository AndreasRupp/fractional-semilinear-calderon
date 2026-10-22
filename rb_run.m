function R=rb_run(d,outdir,p)
root=fileparts(mfilename('fullpath'));addpath(root);
if nargin<2,outdir=fullfile(root,'output',sprintf('example%d',d));end
if nargin<3,p=2;end
if ~exist(outdir,'dir'),mkdir(outdir);end
started=datestr(now,'yyyy-mm-dd HH:MM:SS');exampleClock=tic;
runtimeName='MATLAB';
if exist('OCTAVE_VERSION','builtin'),runtimeName='GNU Octave';end
cfg=rb_config(d,p);
fprintf('Example %d, p = %d: NEW forward solve, maximal Jacobi index %d, started %s.\n',d,cfg.p,cfg.reference_degree,started);
stageClock=tic;ref=rb_reference(cfg);referenceSeconds=toc(stageClock);
fprintf('  Forward solve: %d Newton updates, %d residual checkpoints; residual %.3g; %.6f s.\n', ...
 ref.raw.newtonUpdates,numel(ref.raw.residualHistory),ref.residual,referenceSeconds);
fprintf('  Computing the full exterior response and observation projection...\n');
stageClock=tic;data=rb_data(ref,cfg);dataSeconds=toc(stageClock);
fprintf('  Full flux norm %.10g; data construction %.2f s.\n',data.flux_norm,dataSeconds);
fprintf('  Assembling the finest radial grid (n = %d), observation matrix and SVD...\n',round(1/cfg.noise_h));
stageClock=tic;
fine=rb_prepare(round(1/cfg.noise_h),ref,cfg,data);
fineSeconds=toc(stageClock);
fprintf('  Finest-grid assembly and factorization: %.2f s.\n',fineSeconds);
fprintf('  Actual sizes: state %d, observations %d, coefficient cells %d; B is %d-by-%d.\n', ...
 fine.n,numel(data.g),fine.nd,size(fine.B,1),size(fine.B,2));
runCells=cell(1,numel(cfg.seeds)*(numel(cfg.noise_levels)+numel(cfg.mesh_n)));
runIndex=0;anchors=cell(numel(cfg.seeds),1);
inverseSeconds=zeros(1,numel(runCells));meshAssemblySeconds=zeros(size(cfg.mesh_n));
meshDiagnostics=cell(size(cfg.mesh_n));fineDiagnostics=matrix_diagnostics(fine);
try
fprintf('  Noise study: computing %d inverse reconstructions.\n',numel(cfg.seeds)*numel(cfg.noise_levels));
for k=1:numel(cfg.noise_levels)
 for j=1:numel(cfg.seeds)
  runClock=tic;
  out=rb_inverse(fine,cfg.noise_levels(k),cfg.seeds(j));
  out.study='noise';out.parameter_rule='discrepancy';
  assert(out.discrepancy_attained,'Noise-study discrepancy was not attained.');
  runIndex=runIndex+1;runCells{runIndex}=out;
  inverseSeconds(runIndex)=toc(runClock);
  report_run(runIndex,numel(runCells),out,inverseSeconds(runIndex));
  if out.delta==cfg.mesh_delta,anchors{j}=out;end
 end
end
fprintf('  Mesh study: assembling %d grids and computing %d inverse reconstructions.\n',numel(cfg.mesh_n),numel(cfg.seeds)*numel(cfg.mesh_n));
for k=1:numel(cfg.mesh_n)
 n=cfg.mesh_n(k);
 fprintf('  Assembling radial grid n = %d, h = %.4g...\n',n,1/n);
 stageClock=tic;model=rb_prepare(n,ref,cfg,data);meshAssemblySeconds(k)=toc(stageClock);
 meshDiagnostics{k}=matrix_diagnostics(model);
 fprintf('  Grid assembly and factorization: %.2f s.\n',meshAssemblySeconds(k));
 for j=1:numel(cfg.seeds)
  runClock=tic;
  anchor=anchors{j};assert(~isempty(anchor),'Mesh delta must be a noise-study level.');
  model.alphas=anchor.alpha;model.tau=anchor.tau;
  model.post=chol(model.M+model.tau*model.K);model.fixed_beta=anchor.beta;
  out=rb_inverse(model,cfg.mesh_delta,cfg.seeds(j));
  out.study='mesh';out.parameter_rule='fixed_finest_mesh';
  runIndex=runIndex+1;runCells{runIndex}=out;
  inverseSeconds(runIndex)=toc(runClock);
  report_run(runIndex,numel(runCells),out,inverseSeconds(runIndex));
 end
end
catch problem
 partial=struct('complete',false,'d',d,'cfg',cfg,'runs',{runCells(1:runIndex)}, ...
  'reference',ref.raw,'data',data,'failure',problem.message,'failed_after_run',runIndex, ...
  'execution',version,'created',datestr(now,30));
 save(fullfile(outdir,'partial_results.mat'),'partial','-v7');
 rethrow(problem);
end
assert(runIndex==numel(runCells),'The expected number of runs was not collected.');
runs=[runCells{:}];
R=struct('d',d,'cfg',cfg,'runs',runs,'curveDeltas',cfg.curveDeltas, ...
 'truthRadius',linspace(0,cfg.D,2001)','reference',ref.raw,'data',data, ...
 'execution',version,'created',datestr(now,30));
R.truthQ=ref.q(R.truthRadius);
R.metadata=struct('space','radial P0 conforming for s<1/2', ...
 'physical_domain','unit ball','error_domain_radius',cfg.D, ...
 'mesh_width','radial cell width, not annulus diameter', ...
 'noise_convention','delta = ||noise||_W / ||full flux||_W', ...
 'mesh_protocol','alpha,tau,beta selected on finest grid at mesh_delta, then fixed across h for each seed', ...
 'reference_method','weighted Jacobi Galerkin distinct from radial P0 inverse; Poisson lift evaluated by quadrature', ...
 'coefficient_constraints','no sign constraint, known radial symmetry; q recovered on D');
R.provenance=struct('fresh_computation',true,'loaded_numerical_results',false, ...
 'started',started,'completed',datestr(now,'yyyy-mm-dd HH:MM:SS'), ...
 'source_folder',root,'entry_point','rb_run','runtime_name',runtimeName, ...
 'description','Reference, response, matrices, factorizations and all reconstructions computed in this invocation');
R.timing=struct('reference_seconds',referenceSeconds,'data_seconds',dataSeconds, ...
 'finest_assembly_seconds',fineSeconds,'mesh_assembly_seconds',meshAssemblySeconds, ...
 'inverse_seconds',inverseSeconds,'computation_seconds',toc(exampleClock));
R.numerics=struct('finest',fineDiagnostics,'mesh',{meshDiagnostics}, ...
 'forward_solves',1,'newton_updates',ref.raw.newtonUpdates, ...
 'state_cholesky_factorizations',1+numel(cfg.mesh_n), ...
 'singular_value_decompositions',1+numel(cfg.mesh_n), ...
 'weak_cholesky_factorizations',1+numel(cfg.mesh_n)+numel(cfg.mesh_n)*numel(cfg.seeds), ...
 'inverse_calls',numel(runs),'alpha_candidates_per_noise_call',numel(cfg.alphas));
fprintf('  Saving the %d newly computed runs...\n',numel(runs));
save(fullfile(outdir,'results.mat'),'R','-v7');
rb_tables(R,outdir);
fprintf('Saved %d freshly computed runs to %s (%.2f s including saving).\n',numel(runs),outdir,toc(exampleClock));
type(fullfile(outdir,'summary.txt'));
end
function report_run(index,count,out,elapsed)
fprintf('    [%02d/%02d] %-5s h=%.4g delta=%.0e seed=%d alpha=%.3g e_q=%.6f (%.6f s)\n', ...
 index,count,out.study,out.h,out.delta,out.seed,out.alpha,out.relative_q_error,elapsed);
end
function info=matrix_diagnostics(model)
sigma=model.sigma;threshold=max(size(model.B))*eps(max(sigma));
info=struct('n',model.n,'state_unknowns',model.n,'coefficient_cells',model.nd, ...
 'observation_unknowns',size(model.B,1),'singular_values',sigma, ...
 'numerical_rank_threshold',threshold,'numerical_rank',sum(sigma>threshold), ...
 'alpha_grid',model.alphas,'filter_active_modes',sum(sigma>sqrt(model.alphas),1), ...
 'effective_dimension',sum(sigma.^2./(sigma.^2+model.alphas),1), ...
 'state_factor_residual',norm(model.R'*model.R-model.K-model.M,'fro')/norm(model.K+model.M,'fro'), ...
 'weak_factor_residual',norm(model.post'*model.post-model.M-model.tau*model.K,'fro')/norm(model.M+model.tau*model.K,'fro'));
end
