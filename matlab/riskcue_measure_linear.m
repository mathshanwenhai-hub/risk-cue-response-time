function D = riskcue_measure_linear(R)
%RISK CUE Estimate early-time modal growth and frequency from saved history.
% This is a diagnostic regression, not a proof of convergence or periodicity.
H=R.history;
D=struct('growth',NaN,'omega',NaN,'growth_residual',NaN, ...
 'phase_residual',NaN,'nfit',0,'fit_interval',[NaN NaN]);
if height(H)<8,return;end
z1=H.modal1_real+1i*H.modal1_imag;
z2=H.modal2_real+1i*H.modal2_imag;
amp=sqrt(abs(z1).^2+abs(z2).^2);
tcap=min(R.time_final,max(5,min(30,.30*R.num.Tend)));
id=H.time>0 & H.time<=tcap & amp>1e-12 & H.AP<.05 & H.AN<.05;
if sum(id)<8
 k=min(height(H),max(8,round(.25*height(H))));
 id=false(height(H),1);id(2:k)=true;id=id & amp>1e-12;
end
if sum(id)<5,return;end
t=H.time(id);aa=amp(id);
pg=polyfit(t,log(aa),1);
D.growth=pg(1);
D.growth_residual=sqrt(mean((log(aa)-polyval(pg,t)).^2));
if mean(abs(z2(id)).^2)>mean(abs(z1(id)).^2),z=z2(id);else,z=z1(id);end
ph=unwrap(angle(z));pp=polyfit(t,ph,1);
D.omega=abs(pp(1));
D.phase_residual=sqrt(mean((ph-polyval(pp,t)).^2));
D.nfit=numel(t);D.fit_interval=[t(1) t(end)];
end
