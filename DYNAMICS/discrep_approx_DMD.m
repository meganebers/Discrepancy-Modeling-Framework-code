function out = discrep_approx_DMD(t,y,DMD,nOrder,dt)

b = DMD.w\D_H(:,end);
R = DMD.w*diag(b)*exp(DMD.e*dt.');

out = real(R(1:nOrder,:));