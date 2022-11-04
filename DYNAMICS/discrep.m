function out = discrep(t,y,g,addState)

if addState == 1
    tmp = ['@(y)[',g,';0;0]'];
elseif addState == 2
    tmp = ['@(y)[0;',g,';0]'];
else
    tmp = ['@(y)[0;0;',g,']'];
end

discrepancy = str2func(tmp);

dy_d = @(y) discrepancy(y);

out = dy_d(y);