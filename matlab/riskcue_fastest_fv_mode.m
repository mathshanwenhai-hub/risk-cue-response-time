function M = riskcue_fastest_fv_mode(p,num,excludeZero)
%RISK CUE Fastest eigenmode of the actual cell-centered FV Neumann Laplacian.
% Enumerates every discrete cosine mode resolved by the grid.
if nargin<3,excludeZero=false;end
e=riskcue_equilibrium(p);
best=-Inf;bestOmega=NaN;bestM=0;bestN=0;bestLambda=0;
dx=p.Lx/num.Nx;dy=p.Ly/num.Ny;
for m=0:num.Nx-1
 for n=0:num.Ny-1
  if excludeZero && m==0 && n==0,continue;end
  lh=4/dx^2*sin(m*pi/(2*num.Nx))^2+4/dy^2*sin(n*pi/(2*num.Ny))^2;
  A=riskcue_mode(lh,p.tau,p,e);
  z=eig(A);[g,j]=max(real(z));
  if g>best
   best=g;bestOmega=abs(imag(z(j)));bestM=m;bestN=n;bestLambda=lh;
  end
 end
end
M=struct('m',bestM,'n',bestN,'lambda_h',bestLambda, ...
 'growth',best,'omega',bestOmega,'tau',p.tau);
end
