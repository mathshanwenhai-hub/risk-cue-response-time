function e = riskcue_equilibrium(p)
%RISK CUE Positive equilibrium; select the smallest admissible positive root.
% Root scanning also detects multiple positive equilibria instead of assuming
% uniqueness for every parameter set. The benchmark satisfies a*h*K <= 1.
K=p.r/p.r1; F=@(n)p.a*n./(1+p.a*p.h*n);
if p.b*F(K)<=p.delta
 error('riskcue:equilibrium','Predator invasion inequality b*F(K)>delta fails.');
end
Nmin=p.delta/(p.a*(p.b-p.delta*p.h));
x=linspace(Nmin+max(1e-12,1e-9*Nmin),K,1001);
PofN=@(n)(p.b*F(n)-p.delta)/p.delta1;
g=@(n)p.r-p.r1*n-p.a*PofN(n)./(1+p.a*p.h*n);
y=g(x); rootsN=[];
for j=1:numel(x)-1
 if y(j)*y(j+1)<0
  rootsN(end+1)=fzero(g,[x(j),x(j+1)]); %#ok<AGROW>
 end
end
if isempty(rootsN),error('riskcue:equilibrium','No positive coexistence root was bracketed.');end
if numel(rootsN)>1
 warning('riskcue:multipleEquilibria','Multiple positive equilibria; choosing smallest N.');
end
e.N=rootsN(1); e.P=PofN(e.N); e.W=p.gamma*e.P/p.mu;
den=1+p.a*p.h*e.N;
e.f=F(e.N);e.fp=p.a/den^2;
e.fpp=-2*p.a^2*p.h/den^3;e.fppp=6*p.a^3*p.h^2/den^4;
e.beta0=-p.r+2*p.r1*e.N+e.P*e.fp;
e.all_positive_N_roots=rootsN;
end
