function study = run_siap_study(task,options)
%RUN_SIAP_STUDY Reproducible entry point. Default is THEORY ONLY.
%   run_siap_study                 % fast algebra, tables, Figure 1
%   run_siap_study('baseline')     % tau=0,1,8; three nonlinear runs
%   run_siap_study('near_entry')   % centered-flux pure/symmetric branches
%   run_siap_study('near_exit')    % tau=6,6.8,8
%   run_siap_study('refine')       % one tau at baseline and dt/2
%   run_siap_study('continuation') % increasing tau, reusing final states
% Optional overrides, e.g.
%   opt.num=struct('Nx',40,'Ny',40,'Tend',20);
%   run_siap_study('baseline',opt); % SMOKE TEST ONLY, not publication data
%   opt.taus=[0 .5 .65 1 2 6 6.8 8];
% No MATLAB PDE toolbox, signal-processing toolbox, or parallel toolbox used.
if nargin<1 || isempty(task),task='theory';end
if nargin<2,options=struct();end
[p,num]=riskcue_parameters();
root=fileparts(fileparts(mfilename('fullpath')));
outdir=fullfile(root,'results');figdir=fullfile(root,'figures');
if ~exist(outdir,'dir'),mkdir(outdir);end
if ~exist(figdir,'dir'),mkdir(figdir);end
if isfield(options,'p'),p=merge_fields(p,options.p);end
if isfield(options,'num'),num=merge_fields(num,options.num);end
fprintf('Task=%s; grid %dx%d, dt=%g; flux=%s.\n',task,num.Nx,num.Ny,num.dt,num.flux);
study=struct('task',task,'p',p,'num',num,'files',{{}});
switch lower(task)
 case 'theory'
  S=riskcue_spectrum(p,num,outdir);study.spectrum=S;
  ki=[S.entry.m,S.entry.n];ko=[S.exit.m,S.exit.n];
  study.entry=riskcue_normalform(p,ki,S.entry.tau_minus,outdir);
  study.exit=riskcue_normalform(p,ko,S.exit.tau_plus,outdir);
  riskcue_plot('theory',S,figdir);
  save(fullfile(outdir,'theory_bundle.mat'),'study');
  fprintf('Theory complete. No nonlinear PDE simulation has been run.\n');
  return
 case 'baseline'
  taus=[0 1 8];if isfield(options,'taus'),taus=options.taus;end
  jobs=make_jobs(taus,p,num,'broadband');
 case 'near_exit'
  taus=[6 6.8 8];if isfield(options,'taus'),taus=options.taus;end
  num.Tend=max(num.Tend,240);if isfield(options,'num') && isfield(options.num,'Tend'),num.Tend=options.num.Tend;end
  jobs=make_jobs(taus,p,num,'broadband');
 case 'near_entry'
  S=riskcue_spectrum(p,num,outdir);[tc,j]=min(S.modes.tau_minus_FV);
  if ~isfinite(tc),error('riskcue:entry','No grid crossing.');end
  mode=[S.modes.m(j),S.modes.n(j)];
  taus=tc+[.02 .04];if isfield(options,'taus'),taus=options.taus;end
  num.flux='centered';num.mode=mode;num.Tend=max(num.Tend,600);
  num.initial_amplitude=.001;num.store_snapshots=true;
  if isfield(options,'num'),num=merge_fields(num,options.num);end
  jobs=[make_jobs(taus,p,num,'pure'),make_jobs(taus,p,num,'symmetric')];
  fprintf('Continuous entry %.9g, selected FV entry %.9g, mode=(%d,%d).\n', ...
   S.entry.tau_minus,tc,mode(1),mode(2));
  fprintf('Near-critical runs are long. T is a user-controlled horizon, not a convergence guarantee.\n');
 case 'refine'
  jobs=make_jobs(p.tau,p,num,'broadband');
  fine=num;fine.dt=num.dt/2;jobs=[jobs,make_jobs(p.tau,p,fine,'broadband')];
 case 'continuation'
  taus=[.5 .7 1 2 4 6 6.8 8];if isfield(options,'taus'),taus=options.taus;end
  jobs=make_jobs(taus,p,num,'broadband');
 otherwise,error('riskcue:task','Unknown task. See help run_siap_study.');
end
runs=cell(1,numel(jobs));initial=[];
for j=1:numel(jobs)
 q=jobs(j).p;nu=jobs(j).num;
 tag=sprintf('%s_tau%.8g_N%d_dt%.8g_T%.8g_%s_%s', ...
  task,q.tau,nu.Nx,nu.dt,nu.Tend,nu.initial_kind,nu.flux);
 tag=[tag sprintf('_seed%.4g',nu.transverse_seed)];
 tag=strrep(tag,'.','p');fname=fullfile(outdir,[tag '.mat']);
 if ~strcmpi(task,'continuation'),initial=[];end
 runs{j}=riskcue_simulate(q,nu,fname,initial);study.files{j}=fname;
 if strcmpi(task,'continuation')
  if ~strcmp(runs{j}.status,'completed')
   warning('riskcue:continuation','Stopping continuation after failed run.');break;
  end
  initial=struct('P',runs{j}.P,'N',runs{j}.N,'W',runs{j}.W);
 end
end
runs=runs(~cellfun(@isempty,runs));study.runs=runs;
save(fullfile(outdir,[lower(task) '_index.mat']),'study','-v7.3');
if strcmpi(task,'baseline'),riskcue_plot('baseline',runs,figdir);end
if strcmpi(task,'near_entry'),riskcue_plot('branches',runs,figdir);end
end
function v=merge_fields(v,u)
f=fieldnames(u);for j=1:numel(f),v.(f{j})=u.(f{j});end
end
function jobs=make_jobs(taus,p,num,kind)
jobs=struct('p',{},'num',{});
for j=1:numel(taus)
 q=p;q.tau=taus(j);nu=num;nu.initial_kind=kind;
 jobs(end+1)=struct('p',q,'num',nu); %#ok<AGROW>
end
end
