function out = lorenz_approx_SINDy_discrep(t,y,ahat,sigma,beta,rho,polyorder,usesine)

% poolData(x_d,n,polyorder,usesine);

dy = @(y)[
sigma*(y(2)-y(1));
y(1)*(rho-y(3))-y(2);
y(1)*y(2)-beta*y(3);
];

tmp = dy(y);

yPool = poolData(y',length(y),polyorder,usesine);
out = (yPool*ahat)' + tmp;