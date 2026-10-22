function rb_tables(R,outdir)
fid=fopen(fullfile(outdir,'runs.csv'),'w');assert(fid>=0);
fprintf(fid,'dimension,study,seed,h,delta,alpha,tau,beta,q_error,u_error,w_error,visibility,residual,noise_norm,discrepancy_target,discrepancy_attained\n');
for k=1:numel(R.runs)
 z=R.runs(k);
 fprintf(fid,'%d,%s,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%d\n', ...
 R.d,z.study,z.seed,z.h,z.delta,z.alpha,z.tau,z.beta,z.relative_q_error,z.relative_u_error,z.relative_w_error,z.visibility,z.residual,z.noise_norm,z.target,z.discrepancy_attained);
end
fclose(fid);
fid=fopen(fullfile(outdir,'summary.txt'),'w');assert(fid>=0);
fprintf(fid,'Physical dimension %d; p = %d; execution: %s\n',R.d,R.cfg.p,R.execution);
if isfield(R,'provenance')
 fprintf(fid,'Fresh computation started %s; completed %s.\n',R.provenance.started,R.provenance.completed);
 fprintf(fid,'No saved numerical results were loaded. Source folder: %s\n',R.provenance.source_folder);
 fprintf(fid,'Measured computation time before saving: %.3f s.\n',R.timing.computation_seconds);
end
fprintf(fid,'Full projected flux norm %.16g; correction/full ratio %.16g\n',R.data.flux_norm,sqrt(sum(R.data.mass.*R.data.mu.^2))/R.data.flux_norm);
fprintf(fid,'Noise study: h = %.8g, five paired seeds\n',R.cfg.noise_h);
fprintf(fid,'delta   mean e_q   min e_q   max e_q   mean e_u   mean e_w\n');
for delta=R.cfg.noise_levels
 idx=strcmp({R.runs.study},'noise') & [R.runs.delta]==delta;z=R.runs(idx);
 fprintf(fid,'%.2e %.8f %.8f %.8f %.8f %.8f\n',delta,mean([z.relative_q_error]),min([z.relative_q_error]),max([z.relative_q_error]),mean([z.relative_u_error]),mean([z.relative_w_error]));
end
fprintf(fid,'Mesh study: delta = %.2e; finest-grid alpha,tau,beta fixed per seed\n',R.cfg.mesh_delta);
fprintf(fid,'h   mean e_q   min e_q   max e_q   mean e_u   mean e_w\n');
for n=R.cfg.mesh_n
 idx=strcmp({R.runs.study},'mesh') & [R.runs.n]==n;z=R.runs(idx);
 fprintf(fid,'%.8g %.8f %.8f %.8f %.8f %.8f\n',1/n,mean([z.relative_q_error]),min([z.relative_q_error]),max([z.relative_q_error]),mean([z.relative_u_error]),mean([z.relative_w_error]));
end
fprintf(fid,'Coarse mesh residuals may exceed the noise discrepancy: this study fixes the regularization to isolate mesh effects.\n');
fclose(fid);
end
