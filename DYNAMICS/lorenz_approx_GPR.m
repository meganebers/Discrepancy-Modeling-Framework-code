function [out] = lorenz_approx_GPR(t,y,sigma,beta,rho,DMdl)

dy = @(y)[
sigma*(y(2)-y(1));
y(1)*(rho-y(3))-y(2);
y(1)*y(2)-beta*y(3);
];

tmp = dy(y);

D_approx = [];

for i = 1:numel(y)
    D_approx(i) = predict(DMdl{i},y.');
end

%pred = D_approx.';
out = tmp + D_approx.';
