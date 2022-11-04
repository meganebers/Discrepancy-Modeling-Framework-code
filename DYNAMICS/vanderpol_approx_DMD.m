function [out] = vanderpol_approx_DMD(t,y,alpha,DMD,dt,T1)

dy = @(y)[y(2); ...
    alpha*(1 - y(1)^2)*y(2) - y(1)];

tmp = dy(y);

% Calculate next step of discrepancy by "pinning" b_n onto the state space
% by recalculating it with Phi each time step and then adding to d_dy
% R ~ w*diag(b)*exp(e*t')
%b = DMD.w\(tmp-Ef1);

T = dt+T1;
Ef_tmp = DMD.Phi*diag(DMD.b)*exp(DMD.lambda*T);

out =  tmp + real(Ef_tmp(1:2));