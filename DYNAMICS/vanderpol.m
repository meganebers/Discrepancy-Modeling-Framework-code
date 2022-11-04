function out = vanderpol(t,y,mu)

dy = @(y)[ y(2); ...
           mu*(1 - y(1)^2)*y(2) - y(1)];

out = dy(y);