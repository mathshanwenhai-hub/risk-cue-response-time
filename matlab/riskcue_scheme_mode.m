function M = riskcue_scheme_mode(mode,p,num)
%RISK CUE Dominant SBDF2/IMEX amplification factor for one FV cosine mode.
% This predicts the exact linear growth/frequency of the time-discrete scheme
% after the Euler startup transient has decayed.
dx=p.Lx/num.Nx;dy=p.Ly/num.Ny;
lh=4/dx^2*sin(mode(1)*pi/(2*num.Nx))^2+4/dy^2*sin(mode(2)*pi/(2*num.Ny))^2;
e=riskcue_equilibrium(p);dt=num.dt;
alpha0=p.delta1*e.P;eta=e.P*(p.xi*lh+p.b*e.fp);
beta0=-p.r+2*p.r1*e.N+e.P*e.fp;
if p.tau==0
 m=p.mu+p.dW*lh;
 E=[-alpha0,eta;-e.f-p.chi*p.gamma*e.N*lh/m,-beta0];
 Aimp=diag([-p.dP*lh,-p.dN*lh]);I=eye(2);
else
 E=[-alpha0,eta,0;-e.f,-beta0,-p.chi*e.N*lh;p.gamma/p.tau,0,0];
 Aimp=diag([-p.dP*lh,-p.dN*lh,-(p.dW*lh+p.mu)/p.tau]);I=eye(3);
end
D=1.5*I-dt*Aimp;
G=[D\(2*I+2*dt*E),D\(-.5*I-dt*E);I,zeros(size(I))];
z=eig(G);[~,j]=max(abs(z));g=z(j);
M=struct('lambda_h',lh,'amplification',g,'growth',log(abs(g))/dt, ...
 'omega',abs(angle(g))/dt);
end
