function [out] = vanderpol_approx_GPR_discrep(t,y,alpha,DMdl)

dy = @(y)[y(2); ...
          alpha*(1 - y(1)^2)*y(2) - y(1)];

tmp = dy(y);

D_approx = [];

for i = 1:numel(y)
    D_approx(i) = predict(DMdl{i},y.');
end

%pred = D_approx.';
out = tmp + D_approx.';
