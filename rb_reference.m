function ref=rb_reference(cfg)
options=struct('amplitude',cfg.A,'radius',cfg.q_radius,'power',cfg.p, ...
    'inputCentre',cfg.input_center,'inputRadius',cfg.input_width);
if isfield(cfg,'q_amplitude'),options.qAmplitude=cfg.q_amplitude;end
n=256;if isfield(cfg,'reference_degree'),n=cfg.reference_degree;end
options.quadratureOrder=max(512,4*n);
raw=sj_reference(cfg.d,cfg.s,n,options);
ref=struct('raw',raw,'degree',n,'residual',raw.residual,'minimumState',raw.minimumState, ...
    'maximumState',raw.maximumState,'converged',raw.converged,'positiveBranch',raw.positiveBranch);
ref.u0=@(r)sj_harmonic(r,raw);
ref.u=@(r)state(r,raw);
ref.q=@(r)raw.qAmplitude*sj_bump(r/raw.radius);
ref.response=@(rho)full_response(rho,raw);
ref.correction_response=@(rho)correction_response(rho,raw);
ref.background_response=@(rho)background_response(rho,raw);
end
function u=state(r,raw)
u=sj_evaluate(r,raw);
end
function g=full_response(rho,raw)
g=sj_flux(rho,raw);
end
function mu=correction_response(rho,raw)
[~,~,mu]=sj_flux(rho,raw);
end
function g0=background_response(rho,raw)
[~,g0]=sj_flux(rho,raw);
end
