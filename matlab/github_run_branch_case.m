function R = github_run_branch_case(tau,kind,L,Nx,dt,Tend,seed)
%GITHUB_RUN_BRANCH_CASE Clean two-mode Hopf validation on a smaller square.
% L=5 is chosen because the first critical family (2,3)/(3,2) is separated
% from the next response-time crossing, unlike the denser L=10 spectrum.
if nargin<2 || isempty(kind),kind='pure';end
if nargin<3 || isempty(L),L=5;end
if nargin<4 || isempty(Nx),Nx=50;end
if nargin<5 || isempty(dt),dt=.005;end
if nargin<6 || isempty(Tend),Tend=2000;end
if nargin<7 || isempty(seed),seed=1e-8;end
[p,num]=riskcue_parameters();p.Lx=L;p.Ly=L;
num.Nx=Nx;num.Ny=Nx;num.dt=dt;num.Tend=Tend;num.flux='centered';
num.output_dt=.1;num.progress_every=max(1,round(50/dt));
S=riskcue_spectrum(p,num,'');
tc=S.entry.tau_minus;mode=[S.entry.m,S.entry.n];
v=sort(unique(S.modes.tau_minus));
second=v(find(v>tc+1e-10,1,'first'));
if isempty(second),second=Inf;end
if ~(tau>tc && tau<second)
 error('riskcue:branchWindow','tau=%g must lie between first crossing %.12g and next crossing %.12g.',tau,tc,second);
end
NF=riskcue_normalform(p,mode,tc,'');
p.tau=tau;num.mode=mode;num.initial_kind=kind;num.initial_amplitude=.001;
num.transverse_seed=seed;num.store_snapshots=true;num.snapshot_tail=min(30,Tend);
root=fileparts(fileparts(mfilename('fullpath')));
outdir=fullfile(root,'results');figdir=fullfile(root,'figures');
if ~exist(outdir,'dir'),mkdir(outdir);end
if ~exist(figdir,'dir'),mkdir(figdir);end
tag=sprintf('branch_L%.4g_%s_tau%.8g_N%d_dt%.8g_T%.8g',L,kind,tau,Nx,dt,Tend);
tag=strrep(tag,'.','p');outfile=fullfile(outdir,[tag '.mat']);
R=riskcue_simulate(p,num,outfile,[]);
H=R.history;z1=H.modal1_real+1i*H.modal1_imag;z2=H.modal2_real+1i*H.modal2_imag;
id=H.time>=max(.75*R.time_final,R.time_final-400);
meanZ1=mean(abs(z1(id)).^2);meanZ2=mean(abs(z2(id)).^2);
mu=tau-tc;
if strcmpi(kind,'pure')
 predR2=max(0,-real(NF.kappa)*mu/real(NF.a));
 predZ1=predR2;predZ2=0;
 predOmega=NF.omega+imag(NF.kappa)*mu+imag(NF.a)*predR2;
 transRate=predR2*max(real(NF.pure_transverse_per_R2));
 transverse=abs(z2)./max(abs(z1),eps);
else
 predR2=max(0,-real(NF.kappa)*mu/real(NF.synchronous));
 predZ1=predR2/2;predZ2=predR2/2;
 predOmega=NF.omega+imag(NF.kappa)*mu+imag(NF.synchronous)*predR2;
 transRate=predR2*max(real(NF.synchronous_transverse_per_R2));
 transverse=abs(z1-z2)./max(abs(z1+z2),eps);
end
lateTrans=mean(transverse(id));finalTrans=transverse(end);
if isfinite(R.period_diagnostic.candidate_period)
 measOmega=2*pi/R.period_diagnostic.candidate_period;
else,measOmega=NaN;end
T=table(tau,tc,second,L,Nx,dt,Tend,{kind},mode(1),mode(2), ...
 real(NF.kappa),imag(NF.kappa),real(NF.a),imag(NF.a),real(NF.synchronous),imag(NF.synchronous), ...
 predZ1,predZ2,meanZ1,meanZ2,predOmega,measOmega,transRate,lateTrans,finalTrans, ...
 R.period_diagnostic.field_return_residual,R.maxCFL,R.mass_errors(1),R.mass_errors(2),R.mass_errors(3), ...
 {R.status},{R.message}, ...
 'VariableNames',{'tau','tau_c','next_crossing','L','Nx','dt','Tend','kind','mode_m','mode_n', ...
 'Re_kappa','Im_kappa','Re_a','Im_a','Re_sync','Im_sync','pred_mean_z1_sq','pred_mean_z2_sq', ...
 'measured_mean_z1_sq','measured_mean_z2_sq','pred_omega_nonlinear','measured_omega_late', ...
 'pred_transverse_rate','mean_transverse_ratio','final_transverse_ratio','field_return_residual', ...
 'maxCFL','mass_residualP','mass_residualN','mass_residualW','status','message'});
writetable(T,fullfile(outdir,['summary_' tag '.csv']));
riskcue_plot('branches',{R},figdir);
save(outfile,'R','NF','S','T','-v7.3');
if ~strcmp(R.status,'completed'),error('riskcue:branchStopped','%s',R.message);end
fprintf('Branch case %s tau=%g: tc=%.9g, next=%.9g, predicted |z1|^2=%.6g, measured=%.6g.\n', ...
 kind,tau,tc,second,predZ1,meanZ1);
end
