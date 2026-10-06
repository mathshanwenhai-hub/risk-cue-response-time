function R = riskcue_simulate(p,num,outfile,initial)
%RISK CUE SBDF2 finite-volume simulation, no toolbox beyond base MATLAB.
% Defaults in riskcue_parameters. This routine runs ONE parameter point.
% Diffusion is implicit; taxis/reaction/cue source are extrapolated explicitly.
% First-order upwind taxis; second-order time stepping after Euler startup.
% Negative values are NEVER clipped. A failed run is flagged and saved.
if nargin<3,outfile='';end
if nargin<4,initial=[];end
validateattributes(p.tau,{'numeric'},{'scalar','real','nonnegative','finite'});
validateattributes(num.Nx,{'numeric'},{'scalar','integer','>=',4});
validateattributes(num.Ny,{'numeric'},{'scalar','integer','>=',4});
validateattributes(num.dt,{'numeric'},{'scalar','positive','finite'});
validateattributes(num.Tend,{'numeric'},{'scalar','positive','finite'});
e=riskcue_equilibrium(p);Nx=num.Nx;Ny=num.Ny;
dx=p.Lx/Nx;dy=p.Ly/Ny;dA=dx*dy;area=p.Lx*p.Ly;
x=((1:Nx)-.5)*dx;y=((1:Ny)-.5)*dy;[X,Y]=meshgrid(x,y);
nSteps=round(num.Tend/num.dt);
if nSteps<2 || abs(nSteps*num.dt-num.Tend)>1e-9*max(1,num.Tend)
 error('riskcue:timeGrid','Tend must be an integer multiple of the fixed dt.');
end
dt=num.dt;every=max(1,round(num.output_dt/dt));
Lx=lap1(Nx,dx);Ly=lap1(Ny,dy);
Lap=kron(speye(Nx),Ly)+kron(Lx,speye(Ny));I=speye(Nx*Ny);
Dqs=decomposition(p.mu*I-p.dW*Lap,'chol');
DP1=decomposition(I-dt*p.dP*Lap,'chol');DN1=decomposition(I-dt*p.dN*Lap,'chol');
DP2=decomposition(1.5*I-dt*p.dP*Lap,'chol');DN2=decomposition(1.5*I-dt*p.dN*Lap,'chol');
if p.tau>0
 DW1=decomposition((p.tau+dt*p.mu)*I-dt*p.dW*Lap,'chol');
 DW2=decomposition((1.5*p.tau+dt*p.mu)*I-dt*p.dW*Lap,'chol');
end
mode=num.mode;
if any(mode>=min(Nx,Ny)),error('riskcue:mode','Monitored modes are unresolved.');end
phi1=2*cos(mode(1)*pi*X/p.Lx).*cos(mode(2)*pi*Y/p.Ly);
phi2=2*cos(mode(2)*pi*X/p.Lx).*cos(mode(1)*pi*Y/p.Ly);
% Normalize also if non-square/nonstandard modes are selected.
phi1=phi1/sqrt(mean(phi1(:).^2));phi2=phi2/sqrt(mean(phi2(:).^2));
[monitorRow,monitorQ]=monitor_eigen(mode,p,num,e);
if isempty(initial)
 switch lower(num.initial_kind)
  case 'broadband'
   zP=cos(pi*X/p.Lx).*cos(2*pi*Y/p.Ly)+.4*cos(4*pi*X/p.Lx).*cos(6*pi*Y/p.Ly) ...
       +.23*cos(6*pi*X/p.Lx).*cos(4*pi*Y/p.Ly)+.2*cos(7*pi*X/p.Lx).*cos(3*pi*Y/p.Ly);
   zN=.8*cos(2*pi*X/p.Lx).*cos(pi*Y/p.Ly)-.3*cos(4*pi*X/p.Lx).*cos(6*pi*Y/p.Ly) ...
       +.21*cos(6*pi*X/p.Lx).*cos(4*pi*Y/p.Ly);
   zP=zP/max(abs(zP(:)));zN=zN/max(abs(zN(:)));
   P=e.P*(1+num.initial_amplitude*zP);N=e.N*(1+num.initial_amplitude*zN);
  case {'pure','symmetric'}
   if strcmpi(num.initial_kind,'pure'),phi=phi1;else,phi=(phi1+phi2)/sqrt(2);end
   P=e.P+2*num.initial_amplitude*real(monitorQ(1))*phi;
   N=e.N+2*num.initial_amplitude*real(monitorQ(2))*phi;
   P=P+num.transverse_seed*phi2;
   N=N+.6*num.transverse_seed*phi2;
  otherwise,error('riskcue:initial','Unknown initial_kind.');
 end
 W=reshape(Dqs\(p.gamma*P(:)),Ny,Nx);
else
 if ~isequal(size(initial.P),[Ny Nx]) || ~isequal(size(initial.N),[Ny Nx])
  error('riskcue:initial','Continuation fields must match the current grid.');
 end
 P=initial.P;N=initial.N;
 if p.tau>0 && isfield(initial,'W'),W=initial.W;
 else,W=reshape(Dqs\(p.gamma*P(:)),Ny,Nx);end
end
if any(~isfinite([P(:);N(:);W(:)])) || min([P(:);N(:);W(:)])<0
 error('riskcue:initial','Initial fields must be finite and nonnegative.');
end
names={'time','overlap','rho','consumption','encounter_proxy','Ptotal','Ntotal', ...
 'AP','AN','Mlag','Mlag_squared','C_mean_abundance','C_heterogeneity','C_association', ...
 'deltaC','deltaC_abundance','decomposition_residual','minP','minN','minW', ...
 'maxP','maxN','maxW','maxCFL_sofar','mass_residualP','mass_residualN','mass_residualW', ...
 'modal1_real','modal1_imag','modal2_real','modal2_imag'};
capacity=ceil(nSteps/every)+3;data=nan(capacity,numel(names));
[gP,gN,reP,reN,rate]=explicit(P,N,W,p,dx,dy,num.flux);
maxCFL=dt*rate;masserr=[0 0 0];clock=tic;
ks=1;data(ks,:)=metrics(0,P,N,W,e,p,num,Dqs,dA,area,phi1,phi2,monitorRow,maxCFL,masserr);
snapshots=zeros(Ny,Nx,3,0);snapshot_times=[];
status='completed';message='';t=0;Pprev=[];Nprev=[];Wprev=[];
gPprev=[];gNprev=[];rePprev=[];reNprev=[];
for step=1:nSteps
 if dt*rate>num.CFL_limit
  status='stopped';message=sprintf('Outflow taxis CFL %.4g exceeds diagnostic limit %.4g at t=%.6g. Restart with smaller fixed dt.',dt*rate,num.CFL_limit,t);break;
 end
 if step==1
  Pnext=reshape(DP1\(P(:)+dt*gP(:)),Ny,Nx);
  Nnext=reshape(DN1\(N(:)+dt*gN(:)),Ny,Nx);
  if p.tau>0
   Wnext=reshape(DW1\(p.tau*W(:)+dt*p.gamma*P(:)),Ny,Nx);
  else,Wnext=reshape(Dqs\(p.gamma*Pnext(:)),Ny,Nx);end
 else
  Pnext=reshape(DP2\(2*P(:)-.5*Pprev(:)+dt*(2*gP(:)-gPprev(:))),Ny,Nx);
  Nnext=reshape(DN2\(2*N(:)-.5*Nprev(:)+dt*(2*gN(:)-gNprev(:))),Ny,Nx);
  if p.tau>0
   Wnext=reshape(DW2\(p.tau*(2*W(:)-.5*Wprev(:))+dt*p.gamma*(2*P(:)-Pprev(:))),Ny,Nx);
  else,Wnext=reshape(Dqs\(p.gamma*Pnext(:)),Ny,Nx);end
 end
 if any(~isfinite([Pnext(:);Nnext(:);Wnext(:)]))
  status='stopped';message=sprintf('Nonfinite candidate at t=%.6g. Not a blow-up diagnosis.',step*dt);break;
 end
 if min([Pnext(:);Nnext(:);Wnext(:)]) < -num.negative_tolerance
  status='stopped';message=sprintf('Negative candidate at t=%.6g; no clipping. Restart with smaller dt.',step*dt);break;
 end
 % Discrete mass balances use the ACTUAL extrapolated reaction and source.
 if step==1
  residualP=abs(sum(Pnext(:)-P(:))*dA/dt-sum(reP(:))*dA);
  residualN=abs(sum(Nnext(:)-N(:))*dA/dt-sum(reN(:))*dA);
  if p.tau>0
   residualW=abs(p.tau*sum(Wnext(:)-W(:))*dA/dt+p.mu*sum(Wnext(:))*dA-p.gamma*sum(P(:))*dA);
  else,residualW=abs(p.mu*sum(Wnext(:))*dA-p.gamma*sum(Pnext(:))*dA);end
 else
  residualP=abs(sum(1.5*Pnext(:)-2*P(:)+.5*Pprev(:))*dA/dt-sum(2*reP(:)-rePprev(:))*dA);
  residualN=abs(sum(1.5*Nnext(:)-2*N(:)+.5*Nprev(:))*dA/dt-sum(2*reN(:)-reNprev(:))*dA);
  if p.tau>0
   residualW=abs(p.tau*sum(1.5*Wnext(:)-2*W(:)+.5*Wprev(:))*dA/dt+ ...
    p.mu*sum(Wnext(:))*dA-p.gamma*sum(2*P(:)-Pprev(:))*dA);
  else,residualW=abs(p.mu*sum(Wnext(:))*dA-p.gamma*sum(Pnext(:))*dA);end
 end
 masserr=max(masserr,[residualP residualN residualW]);
 Pprev=P;Nprev=N;Wprev=W;gPprev=gP;gNprev=gN;rePprev=reP;reNprev=reN;
 P=Pnext;N=Nnext;W=Wnext;t=step*dt;
 [gP,gN,reP,reN,rate]=explicit(P,N,W,p,dx,dy,num.flux);maxCFL=max(maxCFL,dt*rate);
 if mod(step,every)==0 || step==nSteps
  ks=ks+1;data(ks,:)=metrics(t,P,N,W,e,p,num,Dqs,dA,area,phi1,phi2,monitorRow,maxCFL,masserr);
  if num.store_snapshots && t>=max(0,num.Tend-num.snapshot_tail)
   snapshots(:,:,:,end+1)=cat(3,P,N,W);snapshot_times(end+1)=t; %#ok<AGROW>
  end
 end
 if num.progress_every>0 && mod(step,num.progress_every)==0
  fprintf('tau=%g, t=%g/%g, min=%.4g, maxP=%.4g, CFLmax=%.4g, elapsed %.1fs\n', ...
   p.tau,t,num.Tend,min([P(:);N(:);W(:)]),max(P(:)),maxCFL,toc(clock));
 end
end
if data(ks,1)<t
 ks=ks+1;data(ks,:)=metrics(t,P,N,W,e,p,num,Dqs,dA,area,phi1,phi2,monitorRow,maxCFL,masserr);
end
R=struct('p',p,'num',num,'equilibrium',e,'status',status,'message',message, ...
 'time_final',t,'wall_seconds',toc(clock),'x',x,'y',y,'P',P,'N',N,'W',W, ...
 'history',array2table(data(1:ks,:),'VariableNames',names),'maxCFL',maxCFL,'mass_errors',masserr, ...
 'snapshots',snapshots,'snapshot_times',snapshot_times);
% Time averages are provisional until windows/refinement agree.
H=R.history;
if height(H)>2 && t>0 && sum(H.time>=(1-num.late_fraction)*t)>1
 t1=(1-num.late_fraction)*t;id=H.time>=t1;
 av=@(v)trapz(H.time(id),v(id))/(H.time(find(id,1,'last'))-H.time(find(id,1)));
 R.average_interval=[H.time(find(id,1)),t];
 R.averages=struct('O',av(H.overlap),'C',av(H.consumption),'Ntotal',av(H.Ntotal), ...
  'Ptotal',av(H.Ptotal),'AP',av(H.AP),'AN',av(H.AN),'Mlag',av(H.Mlag), ...
  'C_abundance',av(H.C_mean_abundance),'C_heterogeneity',av(H.C_heterogeneity), ...
  'C_association',av(H.C_association), ...
  'modal_energy1',av(H.modal1_real.^2+H.modal1_imag.^2), ...
  'modal_energy2',av(H.modal2_real.^2+H.modal2_imag.^2));
else,R.average_interval=[NaN NaN];R.averages=struct();end
R.period_diagnostic=riskcue_period_diagnostic(R);
if ~isempty(outfile)
 folder=fileparts(outfile);if ~isempty(folder) && ~exist(folder,'dir'),mkdir(folder);end
 save(outfile,'R','-v7.3');[folder,base]=fileparts(outfile);
 writetable(R.history,fullfile(folder,[base '_history.csv']));
end
if strcmp(status,'stopped'),warning('riskcue:stopped','%s',message);end
fprintf('Finished tau=%g: status=%s, t=%g, wall=%.1fs\n',p.tau,status,t,R.wall_seconds);
end

function L=lap1(n,h)
e=ones(n,1);L=spdiags([e,-2*e,e],[-1 0 1],n,n);
L(1,1)=-1;L(end,end)=-1;L=L/h^2;
end
function [gp,gn,rp,rn,rate]=explicit(P,N,W,p,dx,dy,flux)
F=p.a*N./(1+p.a*p.h*N);
[tp,cp]=taxis(P,N,p.xi,dx,dy,flux);[tn,cn]=taxis(N,W,-p.chi,dx,dy,flux);
rp=P.*(p.b*F-p.delta-p.delta1*P);rn=p.r*N-p.r1*N.^2-F.*P;
gp=tp+rp;gn=tn+rn;rate=max(cp,cn);
end
function [rhs,rate]=taxis(U,S,coef,dx,dy,flux)
[ny,nx]=size(U);vx=zeros(ny,nx+1);vy=zeros(ny+1,nx);
vx(:,2:nx)=coef*(S(:,2:end)-S(:,1:end-1))/dx;
vy(2:ny,:)=coef*(S(2:end,:)-S(1:end-1,:))/dy;
fx=zeros(size(vx));fy=zeros(size(vy));
if strcmpi(flux,'centered')
 v=vx(:,2:nx);fx(:,2:nx)=v.*(U(:,1:end-1)+U(:,2:end))/2;
 v=vy(2:ny,:);fy(2:ny,:)=v.*(U(1:end-1,:)+U(2:end,:))/2;
elseif strcmpi(flux,'upwind')
 v=vx(:,2:nx);fx(:,2:nx)=max(v,0).*U(:,1:end-1)+min(v,0).*U(:,2:end);
 v=vy(2:ny,:);fy(2:ny,:)=max(v,0).*U(1:end-1,:)+min(v,0).*U(2:end,:);
else,error('riskcue:flux','Unknown face flux.');end
rhs=-(fx(:,2:end)-fx(:,1:end-1))/dx-(fy(2:end,:)-fy(1:end-1,:))/dy;
out=(max(vx(:,2:end),0)-min(vx(:,1:end-1),0))/dx+ ...
    (max(vy(2:end,:),0)-min(vy(1:end-1,:),0))/dy;
rate=max(out(:));
end
function [l,q]=monitor_eigen(mode,p,num,e)
lh=4/(p.Lx/num.Nx)^2*sin(mode(1)*pi/(2*num.Nx))^2+ ...
   4/(p.Ly/num.Ny)^2*sin(mode(2)*pi/(2*num.Ny))^2;
A=riskcue_mode(lh,p.tau,p,e);[V,D]=eig(A);z=diag(D);
positive=find(imag(z)>1e-9);
if isempty(positive),[~,j]=max(real(z));else,[~,k]=max(real(z(positive)));j=positive(k);end
q=V(:,j);q=q/q(1);[VL,DL]=eig(A.');[~,k]=min(abs(diag(DL)-z(j)));
l=VL(:,k).';l=l/(l*q);
if p.tau==0,q=[q;p.gamma*q(1)/(p.mu+p.dW*lh)];l=[l,0];end
end
function v=metrics(t,P,N,W,e,p,num,Dqs,dA,area,phi1,phi2,l,cfl,err)
pm=mean(P(:));nm=mean(N(:));pc=P-pm;nc=N-nm;
AP=sqrt(mean(pc(:).^2))/max(pm,eps);AN=sqrt(mean(nc(:).^2))/max(nm,eps);
O=sum(P(:).*N(:))/sqrt(max(sum(P(:).^2)*sum(N(:).^2),realmin));
rho=NaN;
if AP>num.correlation_floor && AN>num.correlation_floor
 rho=sum(pc(:).*nc(:))/sqrt(sum(pc(:).^2)*sum(nc(:).^2));
end
F=p.a*N./(1+p.a*p.h*N);Fm=p.a*nm/(1+p.a*p.h*nm);
C=sum(P(:).*F(:))*dA;E=sum(P(:).*N(:))*dA;
Cmean=area*pm*Fm;Cheter=area*pm*(mean(F(:))-Fm);
Cassoc=area*mean(pc(:).*(F(:)-mean(F(:))));Cstar=area*e.P*e.f;
Wqs=reshape(Dqs\(p.gamma*P(:)),size(P));
lag=norm(W(:)-Wqs(:))/max(norm(Wqs(:)),realmin);
u1=[mean((P(:)-e.P).*phi1(:));mean((N(:)-e.N).*phi1(:));mean((W(:)-e.W).*phi1(:))];
u2=[mean((P(:)-e.P).*phi2(:));mean((N(:)-e.N).*phi2(:));mean((W(:)-e.W).*phi2(:))];
z1=l*u1;z2=l*u2;
v=[t,O,rho,C,E,area*pm,area*nm,AP,AN,lag,lag^2,Cmean,Cheter,Cassoc, ...
 C-Cstar,Cmean-Cstar,C-Cmean-Cheter-Cassoc,min(P(:)),min(N(:)),min(W(:)), ...
 max(P(:)),max(N(:)),max(W(:)),cfl,err,real(z1),imag(z1),real(z2),imag(z2)];
end
