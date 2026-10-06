function S = riskcue_spectrum(p,num,outdir)
%RISK CUE Complete continuous-PDE mode-window enumeration with proven cutoff.
if nargin<3,outdir='';end
e=riskcue_equilibrium(p);
if e.beta0<0
 error('riskcue:cutoff','Implemented a priori cutoff requires beta0>=0.');
end
C=p.gamma*p.chi*e.N*e.P*(p.xi+p.b*e.fp);
cut=max(1,C/(p.dW*(p.dP+p.dN)^2));
mmax=ceil(p.Lx*sqrt(cut)/pi);nmax=ceil(p.Ly*sqrt(cut)/pi);
rows=[]; tangencies=[];
for m=0:mmax
 for n=0:nmax
  lambda=(m*pi/p.Lx)^2+(n*pi/p.Ly)^2;
  if lambda>cut || (m==0 && n==0),continue;end
  [~,c]=riskcue_mode(lambda,1,p,e);
  if isfinite(c.tau_minus)
   om1=sqrt(c.q0+c.m*c.q1/c.tau_minus);
   om2=sqrt(c.q0+c.m*c.q1/c.tau_plus);
   lh=4/(p.Lx/num.Nx)^2*sin(m*pi/(2*num.Nx))^2+ ...
      4/(p.Ly/num.Ny)^2*sin(n*pi/(2*num.Ny))^2;
   if m<num.Nx && n<num.Ny
    [~,ch]=riskcue_mode(lh,1,p,e);
    hm=ch.tau_minus;hp=ch.tau_plus;
   else,hm=NaN;hp=NaN;end
   rows(end+1,:)=[m,n,lambda,c.tau_minus,c.tau_plus,om1,om2,lh,hm,hp]; %#ok<AGROW>
  elseif isfinite(c.tau_tangent)
   tangencies(end+1,:)=[m,n,lambda,c.tau_tangent]; %#ok<AGROW>
  end
 end
end
if isempty(rows),error('riskcue:noWindows','No strict windows for these parameters.');end
rows=sortrows(rows,4);
T=array2table(rows,'VariableNames',{'m','n','lambda','tau_minus','tau_plus', ...
 'omega_minus','omega_plus','lambda_FV','tau_minus_FV','tau_plus_FV'});
merged=[];
for j=1:size(rows,1)
 lo=rows(j,4);hi=rows(j,5);
 if isempty(merged) || lo>merged(end,2)
  merged(end+1,:)=[lo hi]; %#ok<AGROW>
 else,merged(end,2)=max(merged(end,2),hi);end
end
[~,ji]=min(T.tau_minus);[~,jo]=max(T.tau_plus);
S=struct('p',p,'e',e,'num',num,'cutoff',cut,'modes',T,'merged',merged, ...
 'tangencies',tangencies,'entry',T(ji,:),'exit',T(jo,:));
% Maximal spectral real part over all domain modes under the proven cutoff,
% including the spatially constant mode. Stable high modes cannot destabilize.
lams=unique([0;T.lambda]);
S.taus=linspace(0,9,241)';S.smax=zeros(size(S.taus));
S.omega=zeros(size(S.taus));S.lambda_star=zeros(size(S.taus));
for j=1:numel(S.taus)
 best=-Inf;om=0;ls=0;
 for l=lams'
  A=riskcue_mode(l,S.taus(j),p,e);z=eig(A);[v,id]=max(real(z));
  if v>best,best=v;om=abs(imag(z(id)));ls=l;end
 end
 S.smax(j)=best;S.omega(j)=om;S.lambda_star(j)=ls;
end
% For smax below zero, the above set excludes never-unstable nonzero modes.
% Add all remaining modes below the cutoff for an honest global spectral max.
alll=[];
for m=0:mmax
 for n=0:nmax
  l=(m*pi/p.Lx)^2+(n*pi/p.Ly)^2;
  if l<=cut,alll(end+1)=l;end %#ok<AGROW>
 end
end
alll=unique(alll);
for j=1:numel(S.taus)
 for l=alll
  A=riskcue_mode(l,S.taus(j),p,e);z=eig(A);[v,id]=max(real(z));
  if v>S.smax(j),S.smax(j)=v;S.omega(j)=abs(imag(z(id)));S.lambda_star(j)=l;end
 end
end
S.tail_test_passed=false(size(S.taus));
for j=1:numel(S.taus)
 S.tail_test_passed(j)=tail_check(S.smax(j),S.taus(j),cut,p,e);
end
if ~all(S.tail_test_passed)
 warning('riskcue:spectralTail','Stable abscissae may be truncated where the tail test fails; inspect tail_test_passed.');
end
if ~isempty(outdir)
 if ~exist(outdir,'dir'),mkdir(outdir);end
 writetable(T,fullfile(outdir,'mode_windows_matlab.csv'));
 writetable(array2table(merged,'VariableNames',{'tau_in','tau_out'}),fullfile(outdir,'merged_windows.csv'));
 writetable(table(S.taus,S.smax,S.omega,S.lambda_star,'VariableNames', ...
  {'tau','spectral_abscissa','omega','lambda_star'}),fullfile(outdir,'tau_scan.csv'));
 save(fullfile(outdir,'spectral_results.mat'),'S');
end
fprintf('Equilibrium P*=%.10f, N*=%.10f; cutoff=%.6g\n',e.P,e.N,cut);
disp(S.entry);disp(S.exit);disp(merged);
end

function ok=tail_check(s,tau,cut,p,e)
% Floating-point check of a rigorous sufficient polynomial positivity test.
% Shift z=y+s. Positive coefficients in lambda-cut imply Hurwitz stability
% for every lambda>=cut, hence no omitted mode exceeds the plotted s.
a0=p.delta1*e.P;b0=e.beta0;et=p.xi*e.P;eb=p.b*e.P*e.fp;
q1=[a0+b0,p.dP+p.dN];
q0=[a0*b0+eb*e.f,p.dP*b0+p.dN*a0+et*e.f,p.dP*p.dN];
m=[p.mu,p.dW];R=[0,p.gamma*p.chi*e.N*eb,p.gamma*p.chi*e.N*et];
if tau==0
 b1=padd(q1,2*s);
 b0=padd(conv(m,padd(padd(q0,s*q1),s^2)),R);
 polys={b1,b0};
else
 a1=padd(q1,m/tau);a2=padd(q0,conv(m,q1)/tau);a3=padd(conv(m,q0),R)/tau;
 b1=padd(a1,3*s);b2=padd(padd(a2,2*s*a1),3*s^2);
 b3=padd(padd(padd(a3,s*a2),s^2*a1),s^3);
 h=padd(conv(b1,b2),-b3);polys={b1,b2,b3,h};
end
ok=true;
for h=1:numel(polys)
 a=polys{h};b=zeros(size(a));
 for k=0:numel(a)-1
  for j=k:numel(a)-1,b(k+1)=b(k+1)+a(j+1)*nchoosek(j,k)*cut^(j-k);end
 end
 if any(b<=1e-12),ok=false;end
end
end
function c=padd(a,b)
c=zeros(1,max(numel(a),numel(b)));c(1:numel(a))=a;c(1:numel(b))=c(1:numel(b))+b;
end
