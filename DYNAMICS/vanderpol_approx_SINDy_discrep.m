function out = vanderpol_approx_SINDy_discrep(t,y,ahat,alpha,polyorder,usesine)

% poolData(x_d,n,polyorder,usesine);

dy = @(y)[y(2); ...
          alpha*(1 - y(1)^2)*y(2) - y(1)];

tmp = dy(y);

yPool = poolData(y',length(y),polyorder,usesine);
out = (yPool*ahat)' + tmp;