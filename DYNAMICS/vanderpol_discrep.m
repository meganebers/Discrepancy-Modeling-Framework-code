function out = vanderpol_discrep(t,y,alpha,g)

%g = ['-',num2str(ep),'*y(1).*y(1).*y(1)']; % epsilon discrepancy 

tmp = ['@(y)[0;',g,']'];

discrepancy = str2func(tmp);

vanderpol = @(y)[ y(2); ...
                  alpha*(1 - y(1)^2)*y(2) - y(1)];

dy_g = @(y) vanderpol(y) + discrepancy(y);

out = dy_g(y);