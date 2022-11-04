function [out] = burgers_approx_DMD(t,u0t,k,eps,DMD,dt,T1)

u0=ifft(u0t); u0x=ifft(i*k.*u0t);
rhs = -eps*(k.^2).*u0t - fft(u0x.*u0);

T = dt+T1;
Ef_tmp = DMD.Phi*diag(DMD.b)*exp(DMD.lambda*T);

out = fft(Ef_tmp) + rhs;