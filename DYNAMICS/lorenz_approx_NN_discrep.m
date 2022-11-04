function out = lorenz_approx_NN_discrep(t,y,sigma,beta,rho,net)

lorenz = @(y)[sigma*(y(2)-y(1));
              y(1)*(rho-y(3))-y(2);
              y(1)*y(2)-beta*y(3)];

tmp = lorenz(y);

D_approx = net(y);

out = tmp + D_approx;