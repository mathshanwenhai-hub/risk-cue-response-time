function D = riskcue_period_diagnostic(R)
%RISK CUE Candidate period from modal phase plus a full-field return residual.
% This diagnostic neither proves periodicity nor classifies chaos.
D=struct('status','insufficient_data','candidate_period',NaN, ...
 'phase_fit_residual',NaN,'relative_amplitude_drift',NaN,'field_return_residual',NaN);
H=R.history;if height(H)<20,return;end
id=H.time>=max(.75*R.time_final,R.time_final-20);t=H.time(id);
z=H.modal1_real(id)+1i*H.modal1_imag(id);
z2=H.modal2_real(id)+1i*H.modal2_imag(id);
if mean(abs(z2).^2)>mean(abs(z).^2),z=z2;end
if numel(t)<10 || sqrt(mean(abs(z).^2))<1e-8,return;end
ph=unwrap(angle(z));fit=polyfit(t,ph,1);
if abs(fit(1))<1e-5,return;end
T=2*pi/abs(fit(1));D.candidate_period=T;
D.phase_fit_residual=sqrt(mean((ph-polyval(fit,t)).^2));
nh=floor(numel(z)/2);
D.relative_amplitude_drift=abs(mean(abs(z(1:nh)))-mean(abs(z(nh+1:end))))/max(mean(abs(z)),eps);
D.status='modal_only';
if isempty(R.snapshots) || numel(R.snapshot_times)<10,return;end
st=R.snapshot_times(:);valid=st+T<=st(end);
if sum(valid)<5 || st(end)-st(1)<3*T,return;end
% Rows are time, columns are all three field components.
V=reshape(R.snapshots,[],numel(st)).';
shifted=interp1(st,V,st(valid)+T,'linear');base=V(valid,:);
e=R.equilibrium;ny=R.num.Ny;nx=R.num.Nx;
eq=[repmat(e.P,1,ny*nx),repmat(e.N,1,ny*nx),repmat(e.W,1,ny*nx)];
pert=bsxfun(@minus,base,eq);
D.field_return_residual=norm(shifted-base,'fro')/max(norm(pert,'fro'),eps);
D.status='field_and_modal_candidate';
end
