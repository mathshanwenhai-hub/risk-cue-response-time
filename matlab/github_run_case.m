function R = github_run_case(tau,kind,Nx,dt,Tend,flux,seed)
%GITHUB_RUN_CASE One reproducible CI parameter point with CSV summary.
if nargin<2 || isempty(kind),kind='window';end

function M=mode_prediction(mode,p,num)
e=riskcue_equilibrium(p);
dx=p.Lx/num.Nx;dy=p.Ly/num.Ny;
lh=4/dx^2*sin(mode(1)*pi/(2*num.Nx))^2+4/dy^2*sin(mode(2)*pi/(2*num.Ny))^2;
A=riskcue_mode(lh,p.tau,p,e);z=eig(A);[g,j]=max(real(z));
M=struct('lambda_h',lh,'growth',g,'omega',abs(imag(z(j))));
end
if nargin<3 || isempty(Nx),Nx=100;end
if nargin<4 || isempty(dt),dt=.005;end
if nargin<5 || isempty(Tend),Tend=180;end
if nargin<6 || isempty(flux),flux='upwind';end
if nargin<7 || isempty(seed),seed=0;end
[p,num]=riskcue_parameters();p.tau=tau;
num.Nx=Nx;num.Ny=Nx;num.dt=dt;num.Tend=Tend;num.flux=flux;
num.progress_every=max(0,round(20/dt));
num.output_dt=max(.05,min(.1,20*dt));
Mglobal=riskcue_fastest_fv_mode(p,num,false);
switch lower(kind)
 case 'window'
  Mtarget=riskcue_fastest_fv_mode(p,num,true);
  num.initial_kind='pure';num.mode=[Mtarget.m Mtarget.n];num.initial_amplitude=.001;
 case {'ecology','smoke'}
  num.initial_kind='broadband';num.mode=[4 6];
 case {'pure','symmetric'}
  num.initial_kind=kind;num.mode=[4 6];num.initial_amplitude=.001;
  num.transverse_seed=seed;num.store_snapshots=true;num.snapshot_tail=min(40,Tend);
 otherwise
  error('riskcue:ciKind','Unknown CI case kind.');
end
Pred=mode_prediction(num.mode,p,num);
root=fileparts(fileparts(mfilename('fullpath')));
outdir=fullfile(root,'results');figdir=fullfile(root,'figures');
if ~exist(outdir,'dir'),mkdir(outdir);end
if ~exist(figdir,'dir'),mkdir(figdir);end
tag=sprintf('%s_tau%.8g_N%d_dt%.8g_T%.8g_%s',kind,tau,Nx,dt,Tend,flux);
tag=strrep(tag,'.','p');
outfile=fullfile(outdir,[tag '.mat']);
R=riskcue_simulate(p,num,outfile,[]);
D=riskcue_measure_linear(R);
H=R.history;last=height(H);
if isempty(fieldnames(R.averages))
 avgC=NaN;avgN=NaN;avgO=NaN;avgLag=NaN;ca=NaN;ch=NaN;cc=NaN;
else
 avgC=R.averages.C;avgN=R.averages.Ntotal;avgO=R.averages.O;avgLag=R.averages.Mlag;
 ca=R.averages.C_abundance;ch=R.averages.C_heterogeneity;cc=R.averages.C_association;
end
T=table(tau,Nx,dt,Tend,{kind},{flux},Mglobal.m,Mglobal.n,Mglobal.lambda_h,Mglobal.growth,Mglobal.omega, ...
 num.mode(1),num.mode(2),Pred.lambda_h,Pred.growth,Pred.omega, ...
 D.growth,D.omega,D.growth-Pred.growth,D.omega-Pred.omega,D.nfit, ...
 R.time_final,R.wall_seconds,R.maxCFL,R.mass_errors(1),R.mass_errors(2),R.mass_errors(3), ...
 H.minP(last),H.minN(last),H.minW(last),avgC,avgN,avgO,avgLag,ca,ch,cc,{R.status},{R.message}, ...
 'VariableNames',{'tau','Nx','dt','Tend','kind','flux','global_mode_m','global_mode_n','global_lambda_h', ...
 'global_growth','global_omega','mode_m','mode_n','lambda_h','pred_growth','pred_omega', ...
 'measured_growth','measured_omega','growth_error','omega_error','fit_points','time_final', ...
 'wall_seconds','maxCFL','mass_residualP','mass_residualN','mass_residualW','minP','minN','minW', ...
 'avg_consumption','avg_Ntotal','avg_overlap','avg_cue_mismatch','avg_C_abundance', ...
 'avg_C_heterogeneity','avg_C_association','status','message'});
writetable(T,fullfile(outdir,['summary_' tag '.csv']));
if strcmpi(kind,'smoke')
 % Keep CI smoke tests minimal: no long-time averages or publication figures.
elseif any(strcmpi(kind,{'window','ecology'}))
 riskcue_plot('baseline',{R},figdir);
else
 riskcue_plot('branches',{R},figdir);
end
if ~strcmp(R.status,'completed')
 error('riskcue:ciStopped','Simulation stopped: %s',R.message);
end
end
