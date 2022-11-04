function out = vanderpol_approx_NN_discrep(t,y,alpha,net)

dy = @(y)[y(2); ...
          alpha*(1 - y(1)^2)*y(2) - y(1)];

tmp = dy(y);

D_approx = net(y);

out = tmp + D_approx;