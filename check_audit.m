function audit=check_audit
root=fileparts(mfilename('fullpath'));addpath(root);
folder=fullfile(root,'verification');if ~exist(folder,'dir'),mkdir(folder);end
runtimeName='MATLAB';if exist('OCTAVE_VERSION','builtin'),runtimeName='GNU Octave';end
fprintf('Fresh audit computation in %s %s. No result files are loaded.\n',runtimeName,version);
rows=cell(1,2);clock=tic;
for d=1:2
 cfg=rb_config(d);ref=rb_reference(cfg);data=rb_data(ref,cfg);
 model=rb_prepare(round(1/cfg.noise_h),ref,cfg,data);
 delta=min(cfg.noise_levels);seed=cfg.seeds(1);baseline=rb_inverse(model,delta,seed);
 diagnostic=model;diagnostic.q=2*model.q+1;diagnostic.u=2*model.u+1;diagnostic.wtrue=2*model.wtrue-1;
 mutatedTruth=rb_inverse(diagnostic,delta,seed);
 inverseFields={'q','w','correction','alpha','tau','beta','visibility','noise','residual'};
 for j=1:numel(inverseFields)
  key=inverseFields{j};assert(isequal(baseline.(key),mutatedTruth.(key)), ...
   'Diagnostic truth mutation changed inverse field %s.',key);
 end
 assert(baseline.relative_q_error~=mutatedTruth.relative_q_error,'Diagnostic error did not respond to changed truth.');
 measured=model;direction=model.B*ones(model.n,1);direction=direction/sqrt(sum(data.mass.*direction.^2));
 measured.data.mu=data.mu+.01*data.flux_norm*direction;
 mutatedData=rb_inverse(measured,delta,seed);
 dataChange=norm(mutatedData.q-baseline.q);
 assert(dataChange>64*eps(max(1,norm(baseline.q))),'Changed measurements did not change the reconstruction.');
 physicalCfg=cfg;physicalCfg.q_amplitude=1.2*cfg.q_amplitude;
 physicalRef=rb_reference(physicalCfg);physicalData=rb_data(physicalRef,physicalCfg);
 physicalModel=rb_prepare(model.n,physicalRef,physicalCfg,physicalData);
 mutatedPhysical=rb_inverse(physicalModel,delta,seed);
 physicalFluxChange=sqrt(sum(data.mass.*(physicalData.g-data.g).^2));
 physicalQChange=norm(mutatedPhysical.q-baseline.q);
 assert(physicalFluxChange>64*eps(max(1,data.flux_norm)),'Changed physical coefficient did not change forward data.');
 assert(physicalQChange>64*eps(max(1,norm(baseline.q))),'Changed physical coefficient did not change reconstruction.');
 alphas=unique([baseline.alpha,10.^(-12:-1:-18)]);sensitivity=cell(size(alphas));
 for j=1:numel(alphas)
  probe=model;probe.alphas=alphas(j);value=rb_inverse(probe,delta,seed);
  y=sqrt(data.mass).*(data.mu+value.noise);z=model.Us'*y;
  predicted=sqrt(sum((alphas(j)*z./(model.sigma.^2+alphas(j))).^2)+norm(y-model.Us*z)^2);
  direct=sqrt(sum(data.mass.*(model.B*value.correction-data.mu-value.noise).^2));
  normal=model.B'*(data.mass.*(model.B*value.correction-data.mu-value.noise))+alphas(j)*(model.K+model.M)*value.correction;
  post=(model.M+value.tau*model.K)*value.w-model.K*value.correction;
  sensitivity{j}=struct('alpha',alphas(j),'filter_scale',sqrt(alphas(j)), ...
   'active_modes',sum(model.sigma>sqrt(alphas(j))), ...
   'effective_dimension',sum(model.sigma.^2./(model.sigma.^2+alphas(j))), ...
   'predicted_residual',predicted,'direct_residual',direct, ...
   'residual_difference',abs(predicted-direct),'normal_residual_absolute',norm(normal), ...
   'weak_residual_relative',norm(post)/max(norm(model.K*value.correction),realmin), ...
   'relative_q_error',value.relative_q_error,'q_change_from_selected',norm(value.q-baseline.q), ...
   'q',value.q,'correction',value.correction);
 end
 sigma=model.sigma;threshold=max(size(model.B))*eps(max(sigma));
 S=model.K+model.M;F=model.M+model.tau*model.K;
 rawU=model.u0+baseline.correction(1:model.nd);
 reaction=zeros(model.n,1);reaction(1:model.nd)=baseline.q.*sum(model.weights.*rawU.^model.p,2);
 pdeResidual=model.K*baseline.correction+reaction;
 restrictedDual=norm(model.R'\pdeResidual);
 shellRieszInside=sqrt(sum(pdeResidual(1:model.nd).^2./model.mass(1:model.nd)));
 shellRieszOutside=sqrt(sum(pdeResidual(model.nd+1:end).^2./model.mass(model.nd+1:end)));
 row=struct('dimension',d,'runtime',runtimeName,'version',version, ...
  'diagnostic_truth_mutation_identical',true,'measurement_mutation_q_difference',dataChange, ...
  'physical_mutation_flux_difference',physicalFluxChange,'physical_mutation_q_difference',physicalQChange, ...
  'singular_values',sigma,'numerical_rank_threshold',threshold,'numerical_rank',sum(sigma>threshold), ...
  'singular_value_minimum',min(sigma),'singular_value_maximum',max(sigma), ...
  'state_factor_relative_residual',norm(model.R'*model.R-S,'fro')/norm(S,'fro'), ...
  'weak_factor_relative_residual',norm(model.post'*model.post-F,'fro')/norm(F,'fro'), ...
  'stiffness_symmetry_relative_residual',norm(model.K-model.K','fro')/norm(model.K,'fro'), ...
  'mass_symmetry_relative_residual',norm(model.M-model.M','fro')/norm(model.M,'fro'), ...
  'restricted_pde_dual_norm',restrictedDual,'shell_riesz_inside_D',shellRieszInside, ...
  'shell_riesz_outside_D',shellRieszOutside, ...
  'pde_residual_scope','P0 test-space residual with omitted numerical lift load; not a continuum residual bound', ...
  'baseline',baseline,'mutated_data',mutatedData,'mutated_physical',mutatedPhysical, ...
  'diagnostic_q_error',mutatedTruth.relative_q_error,'sensitivity',{sensitivity}, ...
  'forward_newton_updates',ref.raw.newtonUpdates,'physical_forward_newton_updates',physicalRef.raw.newtonUpdates);
 rows{d}=row;save(fullfile(folder,sprintf('audit_d%d.mat',d)),'row','-v7');
 fprintf('d=%d: truth diagnostic mutation leaves inverse arrays identical; measurement change %.6g, physical change %.6g.\n',d,dataChange,physicalQChange);
 fprintf('  All %d singular values saved; numerical rank %d at threshold %.6g.\n',numel(sigma),row.numerical_rank,threshold);
 fprintf('  Selected alpha %.6g; alpha sensitivity probes include values below the production floor.\n',baseline.alpha);
end
audit=struct('runtime',runtimeName,'version',version,'created',datestr(now,30),'dimensions',{rows},'seconds',toc(clock));
save(fullfile(folder,'fresh_audit.mat'),'audit','-v7');
fprintf('Audit saved in %.6f s. These numerical checks are not continuum error bounds.\n',audit.seconds);
end
