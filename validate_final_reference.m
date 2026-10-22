function validate_final_reference
root=fileparts(mfilename('fullpath'));addpath(root);
outdir=fullfile(root,'verification');if ~exist(outdir,'dir'),mkdir(outdir);end
logfile=fullfile(outdir,'final_reference_checks.txt');
if exist(logfile,'file'),delete(logfile);end
diary(logfile);cleanup=onCleanup(@()diary('off'));
fprintf('Runtime %s. Independent degree and quadrature consistency checks.\n',version);
fprintf('Algebraic and sampled strong residuals are not rigorous continuum error bounds.\n');
checkCells=cell(2,1);
for d=1:2
 cfg=rb_config(d);n=cfg.reference_degree;qo=max(512,4*n);
 options=struct('amplitude',cfg.A,'radius',cfg.q_radius,'qAmplitude',cfg.q_amplitude, ...
  'power',cfg.p,'inputCentre',cfg.input_center,'inputRadius',cfg.input_width, ...
  'quadratureOrder',qo,'inputOrder',128);
 ref=rb_reference(cfg);raw=ref.raw;data=rb_data(ref,cfg);W=data.mass;wn=@(v)sqrt(sum(W.*v.^2));
 coarse=sj_reference(d,cfg.s,n-64,options);
 higher=sj_reference(d,cfg.s,n+64,options);
 nonlinearOptions=options;nonlinearOptions.quadratureOrder=max(1536,6*n);
 nonlinear=sj_reference(d,cfg.s,n,nonlinearOptions);nonlinear=flux_rule(nonlinear,qo);
 inputOptions=options;inputOptions.inputOrder=192;inputRefined=sj_reference(d,cfg.s,n,inputOptions);
 fluxRefined=flux_rule(raw,max(1536,6*n));
 names={'degreeMinus64','degreePlus64','nonlinearQuadrature','inputQuadrature','fluxQuadrature'};
 variants={coarse,higher,nonlinear,inputRefined,fluxRefined};
 r=linspace(0,1,4001)';[u,~,q,w]=sj_evaluate(r,raw);
 row=struct('d',d,'fullFluxNorm',data.flux_norm,'correctionFluxNorm',wn(data.mu), ...
  'forwardResidual',ref.residual,'minimumNoiseAbsolute',min(cfg.noise_levels)*data.flux_norm);
 for j=1:numel(variants)
  trial=variants{j};[g,g0,mu]=sj_flux(data.rho,trial);g=data.projection*g;g0=data.projection*g0;mu=data.projection*mu;
  ut=sj_evaluate(r,trial);label=names{j};
  row.([label,'FullFluxDifference'])=wn(data.g-g);
  row.([label,'CorrectionDifference'])=wn(data.mu-mu);
  row.([label,'BackgroundDifference'])=wn(data.g0-g0);
  row.([label,'StateSupDifference'])=max(abs(u-ut));
 end
 cfgHi=cfg;cfgHi.observation_order=10;dataHi=rb_data(ref,cfgHi);
 row.observationQuadratureFluxDifference=wn(data.g-dataHi.g);
 row.observationQuadratureCorrectionDifference=wn(data.mu-dataHi.mu);
 cfgProjection=cfg;cfgProjection.observation_order=20;dataProjection=rb_data(ref,cfgProjection);
 [gp,~,mup]=sj_flux(dataProjection.rho,raw);nodeWeights=dataProjection.projection'*dataProjection.mass;
 row.unresolvedFullResponseL2=sqrt(sum(nodeWeights.*(gp-repelem(data.g,20)).^2));
 row.unresolvedCorrectionResponseL2=sqrt(sum(nodeWeights.*(mup-repelem(data.mu,20)).^2));
 [rr,ww]=rb_gauss(1800,0,1);sphere=2*pi^(d/2)/gamma(d/2);ww=sphere*ww.*rr.^(d-1);
 [P,lambda]=sj_basis(rr,n,d,cfg.s);uu=sj_evaluate(rr,raw);qq=cfg.q_amplitude*sj_bump(rr/cfg.q_radius);
 differentiated=P*(lambda.*raw.coefficients);strong=differentiated+qq.*uu.^cfg.p;
 row.sampledFullDomainStrongResidualL2=sqrt(sum(ww.*strong.^2));
 row.sampledFullDomainStrongResidualSup=max(abs(strong));outside=rr>cfg.D;
 row.sampledStrongResidualOutsideDL2=sqrt(sum(ww(outside).*strong(outside).^2));
 forcing128=rb_kernel(rr,raw.inputNodes,d,cfg.s)*(raw.inputWeights.*raw.inputValues.*raw.inputNodes.^(d-1));
 forcing192=rb_kernel(rr,inputRefined.inputNodes,d,cfg.s)*(inputRefined.inputWeights.*inputRefined.inputValues.*inputRefined.inputNodes.^(d-1));
 row.liftForcingQuadratureDifferenceL2=sqrt(sum(ww.*(forcing128-forcing192).^2));
 row.minimumState=min(u);row.maximumState=max(u);
 threshold=row.minimumNoiseAbsolute/100;
 for j=1:numel(names)
  assert(row.([names{j},'FullFluxDifference'])<threshold,'Forward sensitivity exceeds one percent of the smallest prescribed noise.');
 end
 assert(row.observationQuadratureFluxDifference<threshold);
 fields=fieldnames(row);fprintf('\nd=%d, baseline degree=%d, nonlinear/flux order=%d, input order=128.\n',d,n,qo);
 fprintf('Degree comparisons hold both quadrature orders at %d; nonlinear, input and flux rules are then varied separately.\n',qo);
 for k=2:numel(fields),fprintf('%s %.17g\n',fields{k},row.(fields{k}));end
 fprintf('Strong residual uses Jacobi differentiation and fresh full-domain quadrature, including Omega outside D.\n');
 checkCells{d}=row;
 save(fullfile(outdir,sprintf('final_reference_d%d.mat',d)),'raw','coarse','higher','nonlinear','inputRefined','fluxRefined','cfg','row','r','u','q','w','rr','strong','-v7');
end
checks=vertcat(checkCells{:});save(fullfile(outdir,'final_reference_checks.mat'),'checks','-v7');
end
function raw=flux_rule(raw,order)
[x,w]=sj_quadrature(order,raw.s,raw.d/2-1);
raw.integrationNodes=sqrt((x+1)/2);raw.integrationWeights=2^(-raw.s-raw.d/2-1)*w;
raw.integrationPolynomials=sj_basis(raw.integrationNodes,raw.degree,raw.d,raw.s);
end
