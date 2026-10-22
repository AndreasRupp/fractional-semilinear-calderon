function run_all
root=fileparts(mfilename('fullpath'));addpath(root);
started=datestr(now,'yyyy-mm-dd HH:MM:SS');totalClock=tic;
powers=[1,2,3];
fprintf('\nFRESH COMPUTATION started %s\n',started);
fprintf('Source folder: %s\nMATLAB/Octave: %s\n',root,version);
fprintf('Known symmetry: even in 1D, radial in 2D; h is radial shell width.\n');
fprintf('Each forward reference has 257 coefficients; the finest inverse has 400 state unknowns.\n');
fprintf('Six forward problems and all 330 inverse runs will be computed for p = 1, 2, 3.\n');
fprintf('The reconstruction panels use p = 2; each error panel compares all three powers.\n');
fprintf('Saved results are not used as input. Elapsed times below are measured.\n\n');
results=cell(2,numel(powers));folders=cell(1,2);
for d=1:2
 outdir=fullfile(root,'output',sprintf('example%d',d));
 folders{d}=outdir;
 for k=1:numel(powers)
  powerdir=outdir;
  if powers(k)~=2,powerdir=fullfile(outdir,sprintf('p%d',powers(k)));end
  results{d,k}=rb_run(d,powerdir,powers(k));
 end
 comparison=results(d,:);
 save(fullfile(outdir,'power_comparison.mat'),'comparison','powers','-v7');
end
fprintf('Both numerical studies are complete and saved. Starting figure export.\n');
exported=false(1,2);
for d=1:2
 fprintf('Example %d: exporting figures from the result just computed...\n',d);
 figureClock=tic;
 try
  files=rb_plot(results(d,:),folders{d});
  exported(d)=numel(files)==2&&all(cellfun(@(f)exist(f,'file')==2,files));
  fprintf('Example %d: figure export returned in %.2f s.\n\n',d,toc(figureClock));
 catch problem
  warning('run_all:figureExport','Example %d export failed: %s. Numerical results remain saved in %s.',d,problem.message,folders{d});
 end
end
fprintf('FRESH COMPUTATION complete: %.2f s elapsed.\n',toc(totalClock));
fprintf('Both numerical examples are saved; figure export completed for %d of 2 examples.\n',sum(exported));
fprintf('Results: %s\n',fullfile(root,'output'));
end
