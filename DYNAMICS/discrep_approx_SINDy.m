function out = discrep_approx_SINDy(t,y,ahat,polyorder,usesine)

% poolData(x_d,n,polyorder,usesine);

yPool = poolData(y',length(y),polyorder,usesine);
out = (yPool*ahat)';