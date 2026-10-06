function R = riskcue_normalform(p,mode,tau,outdir)
%RISK CUE Finite-cosine cubic coefficients at an oscillatory crossing.
% No PDE time integration. q_P=1; mean(phi_j^2)=1; physical time.
% Returns a,b,c12 in the resonant two-mode normal form, the normalized
% synchronous coefficient, mean ecological corrections, and consistency tests.
if nargin<4,outdir='';end
if p.Lx~=p.Ly || any(mode<=0) || mode(1)==mode(2)
 error('riskcue:normalform','This routine is for distinct positive swapped modes on a square.');
end
e=riskcue_equilibrium(p);K=3*max(mode);lambda=lam(mode,p);
[A,c]=riskcue_mode(lambda,tau,p,e);
om=sqrt(c.q0+c.m*c.q1/tau);
q=[1;(c.alpha+1i*om)/c.eta;p.gamma/(c.m+1i*tau*om)];
l=[1,c.eta/(c.beta+1i*om), ...
 -tau*p.chi*e.N*lambda/(c.m+1i*tau*om)*c.eta/(c.beta+1i*om)];
l=l/(l*q);
res=max([norm(A*q-1i*om*q),norm(l*A-1i*om*l),abs(l*q-1),abs(l*conj(q))]);
if res>1e-7,error('riskcue:notHopf','Eigenvector residual %g: tau is not the crossing.',res);end
Ap=zeros(3);Ap(3,:)=[-p.gamma/tau^2,0,c.m/tau^2];kap=l*Ap*q;
Q1=zeros(K+1,K+1,3);Q2=Q1;j=fliplr(mode);
Q1(mode(1)+1,mode(2)+1,:)=reshape(2*q,1,1,3);
Q2(j(1)+1,j(2)+1,:)=reshape(2*q,1,1,3);
B=@(U,V)bilinear(U,V,p,e);C=@(U,V,Z)trilinear(U,V,Z,p,e);
hm=@(U,V)resolvent(B(U,conj(V)),0,tau,p,e);
hp=@(U,V)resolvent(B(U,V),2i*om,tau,p,e);
rep=@(U,V).5*C(U,U,conj(V))+B(U,hm(U,V))+.5*B(conj(V),hp(U,U));
project=@(U)l*reshape(U(mode(1)+1,mode(2)+1,:),3,1)/2;
a=project(rep(Q1,Q1));
b=project(C(Q1,Q2,conj(Q2))+B(Q1,hm(Q2,Q2))+ ...
 B(Q2,hm(Q1,Q2))+B(conj(Q2),hp(Q1,Q2)));
c12=project(rep(Q2,Q1));
Qs=(Q1+Q2)/sqrt(2);Ys=rep(Qs,Qs);
cs=l*(reshape(Ys(mode(1)+1,mode(2)+1,:),3,1)+ ...
 reshape(Ys(j(1)+1,j(2)+1,:),3,1))/(2*sqrt(2));
symmetric_error=abs(cs-(a+b+c12)/2);
H=hm(Q1,Q1);h=reshape(H(1,1,:),3,1);
if norm(imag(h))>1e-7,error('riskcue:mean','Mean correction is not real.');end
h=real(h);
KC=e.f*h(1)+e.P*e.fp*h(2)+2*e.fp*real(q(1)*conj(q(2)))+e.P*e.fpp*abs(q(2))^2;
KO=abs(q(1)/e.P-q(2)/e.N)^2;
KL=2/e.W^2*abs(q(3)-p.gamma*q(1)/c.m)^2;
balP=abs(p.b*KC-(p.delta+2*p.delta1*e.P)*h(1)-2*p.delta1*abs(q(1))^2);
balN=abs(KC-(p.r-2*p.r1*e.N)*h(2)+2*p.r1*abs(q(2))^2);
if max([symmetric_error,balP,balN])>1e-7
 error('riskcue:coefficientCheck','A normalization or balance identity failed.');
end
D=b-a;pure=eig([D,c12;conj(c12),conj(D)]);
D=a-b-3*c12;E=a-b+c12;sync=eig([D,E;conj(E),conj(D)]);
R=struct('mode',mode,'tau',tau,'lambda',lambda,'omega',om,'q',q,'left_row',l, ...
 'kappa',kap,'a',a,'b',b,'c12',c12,'synchronous',cs, ...
 'mean_correction',h,'K_C_per_area',KC,'K_N_per_area',h(2), ...
 'K_overlap',KO,'K_lag_squared',KL, ...
 'pure_transverse_per_R2',pure,'synchronous_transverse_per_R2',sync, ...
 'eigen_residual',res,'symmetric_error',symmetric_error,'mass_errors',[balP balN]);
vals=[kap;a;b;c12;cs];names={'kappa';'a';'b';'c12';'synchronous'};
T=table(names,real(vals),imag(vals),'VariableNames',{'coefficient','real','imag'});
disp(T);fprintf('K_C/area = %.12g, K_N/area = %.12g, K_O = %.12g\n',KC,h(2),KO);
if ~isempty(outdir)
 if ~exist(outdir,'dir'),mkdir(outdir);end
 tag=sprintf('m%d_n%d_tau%.6f',mode(1),mode(2),tau);
 writetable(T,fullfile(outdir,['normalform_' tag '.csv']));
 save(fullfile(outdir,['normalform_' tag '.mat']),'R');
end
end

function l=lam(k,p)
l=(k(1)*pi/p.Lx)^2+(k(2)*pi/p.Ly)^2;
end
function [ks,vs]=active(U)
[a,b]=find(any(abs(U)>1e-14,3));ks=[a-1,b-1];vs=zeros(3,numel(a));
for j=1:numel(a),vs(:,j)=reshape(U(a(j),b(j),:),3,1);end
end
function [k,w]=product(a,b)
% Keep duplicates: four entries with weight 1/4 are subsequently accumulated.
x=[abs(a(1)-b(1)),a(1)+b(1)];y=[abs(a(2)-b(2)),a(2)+b(2)];
k=[x(1) y(1);x(1) y(2);x(2) y(1);x(2) y(2)];w=.25*ones(4,1);
end
function U=accumulate(U,k,v)
if any(k>=size(U,1)),error('riskcue:cosineRange','Cosine support exceeded the exact allocation.');end
U(k(1)+1,k(2)+1,:)=U(k(1)+1,k(2)+1,:)+reshape(v,1,1,3);
end
function Y=bilinear(U,V,p,e)
Y=zeros(size(U));[ku,vu]=active(U);[kv,vv]=active(V);
for i=1:size(ku,1)
 u=vu(:,i);la=lam(ku(i,:),p);
 for j=1:size(kv,1)
  v=vv(:,j);lb=lam(kv(j,:),p);s=u(1)*v(2)+v(1)*u(2);
  r=[-2*p.delta1*u(1)*v(1)+p.b*e.fp*s+p.b*e.P*e.fpp*u(2)*v(2); ...
    -e.fp*s-(2*p.r1+e.P*e.fpp)*u(2)*v(2);0];
  [ks,w]=product(ku(i,:),kv(j,:));
  for h=1:4
   lc=lam(ks(h,:),p);d1=(la-lb-lc)/2;d2=(lb-la-lc)/2;
   tx=[-p.xi*(u(1)*v(2)*d1+v(1)*u(2)*d2); ...
        p.chi*(u(2)*v(3)*d1+v(2)*u(3)*d2);0];
   Y=accumulate(Y,ks(h,:),w(h)*(r+tx));
  end
 end
end
end
function Y=trilinear(U,V,Z,p,e)
Y=zeros(size(U));[ku,vu]=active(U);[kv,vv]=active(V);[kz,vz]=active(Z);
for i=1:size(ku,1)
 u=vu(:,i);
 for j=1:size(kv,1)
  v=vv(:,j);
  for h=1:size(kz,1)
   z=vz(:,h);s=u(1)*v(2)*z(2)+v(1)*u(2)*z(2)+z(1)*u(2)*v(2);
   t=u(2)*v(2)*z(2);
   r=[p.b*e.fpp*s+p.b*e.P*e.fppp*t;-e.fpp*s-e.P*e.fppp*t;0];
   [k1,w1]=product(ku(i,:),kv(j,:));
   for h1=1:4
    [k2,w2]=product(k1(h1,:),kz(h,:));
    for h2=1:4,Y=accumulate(Y,k2(h2,:),w1(h1)*w2(h2)*r);end
   end
  end
 end
end
end
function Y=resolvent(U,shift,tau,p,e)
Y=zeros(size(U));[ks,vs]=active(U);
for j=1:size(ks,1)
 A=riskcue_mode(lam(ks(j,:),p),tau,p,e);
 M=shift*eye(3)-A;
 if rcond(M)<1e-12,warning('riskcue:illConditioned','Near-singular homological matrix.');end
 Y=accumulate(Y,ks(j,:),M\vs(:,j));
end
end
