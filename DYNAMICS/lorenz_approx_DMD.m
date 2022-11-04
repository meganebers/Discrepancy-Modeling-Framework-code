function [out] = lorenz_approx_DMD(t,y,sigma,beta,rho,DMD,dt,T1)

dy = @(y)[
sigma*(y(2)-y(1));
y(1)*(rho-y(3))-y(2);
y(1)*y(2)-beta*y(3);
];

tmp = dy(y);

T = dt+T1;
Ef_tmp = DMD.Phi*diag(DMD.b)*exp(DMD.lambda*T);

out =  tmp + real(Ef_tmp(1:3));