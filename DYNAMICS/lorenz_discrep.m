function out = lorenz_discrep(t,y,sigma,beta,rho,g)
%{
if addState == 1
    tmp = ['@(y)[',g,';0;0]'];
elseif addState == 2
    tmp = ['@(y)[0;',g,';0]'];
else
    tmp = ['@(y)[0;0;',g,']'];
end
%}

tmp = ['@(y)[',g,';0;0]'];

discrepancy = str2func(tmp);

lorenz = @(y)[
sigma*(y(2)-y(1));
y(1)*(rho-y(3))-y(2);
y(1)*y(2)-beta*y(3);
];

dy_d = @(y) lorenz(y) + discrepancy(y);

out = dy_d(y);