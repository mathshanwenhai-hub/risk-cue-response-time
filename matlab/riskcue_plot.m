function riskcue_plot(task,data,figdir)
%RISK CUE Export figures only from actual computations passed to this function.
% Uses base MATLAB print('-dpdf','-painters'); no fictitious trajectory data.
if ~exist(figdir,'dir'),mkdir(figdir);end
switch lower(task)
 case 'theory'
  S=data;p=S.p;e=S.e;
  f=figure('Color','w','Position',[100 100 1100 350]);
  subplot(1,3,1);ls=linspace(0,35,350);ts=[0 .5 1 2 8];hold on;
  for t=ts
   ss=zeros(size(ls));
   for j=1:numel(ls),z=eig(riskcue_mode(ls(j),t,p,e));ss(j)=max(real(z));end
   plot(ls,ss,'LineWidth',1.3,'DisplayName',sprintf('\\tau=%g',t));
  end
  plot(ls,zeros(size(ls)),':','HandleVisibility','off');xlabel('\lambda');ylabel('max Re eigenvalue');
  title('Continuous-mode dispersion');legend('Location','best');grid on;
  subplot(1,3,2);hold on;
  T=S.modes;keep=T.m<=T.n;T=T(keep,:);T=sortrows(T,'lambda');
  for j=1:height(T),plot([T.tau_minus(j),T.tau_plus(j)],[T.lambda(j),T.lambda(j)],'LineWidth',1.0);end
  xlabel('\tau');ylabel('\lambda_{mn}');title('Mode-wise strict windows');grid on;
  subplot(1,3,3);plot(S.taus,S.smax,'LineWidth',1.5);hold on;
  plot(S.taus,zeros(size(S.taus)),':');
  yl=ylim;plot([S.entry.tau_minus S.entry.tau_minus],yl,'--');
  plot([S.exit.tau_plus S.exit.tau_plus],yl,'--');
  xlabel('\tau');ylabel('spectral abscissa (tail checked)');title('Entry and exit');grid on;
  export(f,fullfile(figdir,'fig01_windows.pdf'));
 case 'baseline'
  R=data;R=R(cellfun(@(r)strcmp(r.status,'completed'),R));
  if isempty(R),warning('riskcue:noPlot','No completed runs to plot.');return;end
  [~,idx]=sort(cellfun(@(r)r.p.tau,R));R=R(idx);
  if numel(R)>3,Rshow=R(unique(round(linspace(1,numel(R),3))));else,Rshow=R;end
  nr=numel(Rshow);mn=inf(1,3);mx=-inf(1,3);vars={'P','N','W'};
  for j=1:nr
   for v=1:3,U=Rshow{j}.(vars{v});mn(v)=min(mn(v),min(U(:)));mx(v)=max(mx(v),max(U(:)));end
  end
  f=figure('Color','w','Position',[60 60 360*nr 760]);
  for j=1:nr
   for v=1:3
    subplot(3,nr,(v-1)*nr+j);r=Rshow{j};imagesc(r.x,r.y,r.(vars{v}));axis xy equal tight;
    if mx(v)>mn(v),caxis([mn(v),mx(v)]);end
    colorbar;xlabel('x');ylabel('y');title(sprintf('%s, \\tau=%g, t=%g',vars{v},r.p.tau,r.time_final));
   end
  end
  export(f,fullfile(figdir,'fig03_dynamics.pdf'));
  tau=cellfun(@(r)r.p.tau,R);v0=R{1}.p.Lx*R{1}.p.Ly;
  C0=v0*R{1}.equilibrium.P*R{1}.equilibrium.f;N0=v0*R{1}.equilibrium.N;
  C=cellfun(@(r)r.averages.C,R);NC=cellfun(@(r)r.averages.Ntotal,R);
  O=cellfun(@(r)r.averages.O,R);ML=cellfun(@(r)r.averages.Mlag,R);
  da=cellfun(@(r)r.averages.C_abundance-C0,R);
  dh=cellfun(@(r)r.averages.C_heterogeneity,R);dc=cellfun(@(r)r.averages.C_association,R);
  f=figure('Color','w','Position',[80 80 1000 700]);
  subplot(2,2,1);plot(tau,(C-C0)/C0,'o-','LineWidth',1.3);hold on;
  plot(tau,da/C0,'s--',tau,dh/C0,'d--',tau,dc/C0,'^--');grid on;
  xlabel('\tau');ylabel('relative consumption change');legend('total','abundance','heterogeneity','association','Location','best');
  title('Exact decomposition; provisional time means');
  subplot(2,2,2);plot(tau,(NC-N0)/N0,'o-','LineWidth',1.3);grid on;xlabel('\tau');ylabel('relative prey-total change');
  subplot(2,2,3);plot(tau,O,'o-','LineWidth',1.3);grid on;xlabel('\tau');ylabel('mean overlap');
  subplot(2,2,4);plot(tau,ML,'o-','LineWidth',1.3);grid on;xlabel('\tau');ylabel('mean cue mismatch');
  export(f,fullfile(figdir,'fig04_ecology.pdf'));
  f=figure('Color','w','Position',[80 80 1000 650]);
  for j=1:numel(R)
   H=R{j}.history;label=sprintf('\\tau=%g',R{j}.p.tau);
   subplot(2,2,1);hold on;plot(H.time,H.AP,'DisplayName',label);ylabel('A_P');
   subplot(2,2,2);hold on;plot(H.time,H.consumption,'DisplayName',label);ylabel('C');
   subplot(2,2,3);hold on;plot(H.time,H.overlap,'DisplayName',label);ylabel('O');
   subplot(2,2,4);hold on;plot(H.time,H.Mlag,'DisplayName',label);ylabel('M_{lag}');
  end
  for k=1:4,subplot(2,2,k);grid on;xlabel('time');legend('Location','best');end
  export(f,fullfile(figdir,'supp_histories.pdf'));
 case 'branches'
  R=data;R=R(cellfun(@(r)strcmp(r.status,'completed'),R));
  if isempty(R),warning('riskcue:noPlot','No completed branch runs.');return;end
  f=figure('Color','w','Position',[80 80 1000 650]);
  for j=1:numel(R)
   r=R{j};H=r.history;z1=H.modal1_real+1i*H.modal1_imag;z2=H.modal2_real+1i*H.modal2_imag;
   label=sprintf('%s, \\tau=%.5g',r.num.initial_kind,r.p.tau);
   subplot(2,2,1);hold on;plot(H.time,abs(z1),'DisplayName',label);ylabel('|z_1|');
   subplot(2,2,2);hold on;plot(H.time,abs(z2),'DisplayName',label);ylabel('|z_2|');
   subplot(2,2,3);hold on;plot(H.time,H.consumption,'DisplayName',label);ylabel('consumption');
   subplot(2,2,4);hold on;plot(real(z1),imag(z1),'DisplayName',label);xlabel('Re z_1');ylabel('Im z_1');
  end
  for k=1:3,subplot(2,2,k);grid on;xlabel('time');legend('Location','best');end
  subplot(2,2,4);grid on;axis equal;
  export(f,fullfile(figdir,'fig02_branches.pdf'));
 otherwise,error('riskcue:plot','Unknown plot task.');
end
end
function export(f,name)
set(f,'PaperPositionMode','auto');
print(f,name,'-dpdf','-painters');
fprintf('Exported %s\n',name);
end
