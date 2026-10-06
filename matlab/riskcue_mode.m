function [A,c] = riskcue_mode(lambda,tau,p,e)
%RISK CUE Linearization and response-time quadratic for one mode.
validateattributes(lambda,{'numeric'},{'scalar','real','nonnegative','finite'});
validateattributes(tau,{'numeric'},{'scalar','real','nonnegative','finite'});
c.alpha=p.dP*lambda+p.delta1*e.P;
c.beta=p.dN*lambda-p.r+2*p.r1*e.N+e.P*e.fp;
c.eta=e.P*(p.xi*lambda+p.b*e.fp);c.m=p.mu+p.dW*lambda;
c.q1=c.alpha+c.beta;c.q0=c.alpha*c.beta+c.eta*e.f;
c.R=p.gamma*p.chi*e.N*lambda*c.eta;
if c.q1<=0 || c.q0<=0
 error('riskcue:kinetics','The kinetic damping hypothesis fails at lambda=%g.',lambda);
end
c.Rcrit=c.m*c.q1*(c.q1+2*sqrt(c.q0));
c.H=[c.q1*c.q0,c.m*c.q1^2-c.R,c.m^2*c.q1];
c.tau_minus=NaN;c.tau_plus=NaN;c.tau_tangent=NaN;
margin=c.R-c.Rcrit; tol=1e-11*max(1,abs(c.Rcrit));
if margin>tol
 D=c.R-c.m*c.q1^2;
 disc=D^2-4*c.m^2*c.q1^2*c.q0;
 c.tau_plus=(D+sqrt(max(disc,0)))/(2*c.q1*c.q0);
 c.tau_minus=(c.m^2/c.q0)/c.tau_plus;
elseif abs(margin)<=tol
 c.tau_tangent=c.m/sqrt(c.q0);
end
if tau==0
 A=[-c.alpha,c.eta;-e.f-p.chi*p.gamma*e.N*lambda/c.m,-c.beta];
else
 A=[-c.alpha,c.eta,0; -e.f,-c.beta,-p.chi*e.N*lambda; ...
  p.gamma/tau,0,-c.m/tau];
end
end
