function Xi = sparsifyDynamics_SSR(Theta,delta)
% Compute sparse regression of a discrepancy (ONE STATE): iterative least squares
% SINDy with Stepwise Sparse Regressor [1]
% Code written by Megan Auger
% 11/13/19

% [1] Boninsegna, Nuske, Celmenti. Sparse learning of stochastic dynamical equations. J. Chem. Phys. 148, 241723 (2018); https://doi.org/10.1063/1.5018409148

%% Peform sparse regression

% initial guess: Least-squares
Xi = Theta\delta;  
[n,m]=size(Theta);


% get rid of zero terms in Xi and repsectively Theta and dXdt
if any(Xi == 1)
    Xi_hat = find(Xi == 0);
    Theta(:,Xi_hat) = [];    
    Xi = Theta\delta;  % initial guess: Least-squares
    
    M = m - length(Xi_hat);
else
    % if no coeffs with zeros OR once removed, find min(abs(Xi)) = 0 and re-gregress
    M = m;
end

for i = 1:M-1
    Xi(abs(Xi) == min(abs(Xi(:)))) = 0;
    %[row, col] = find(Xi == 0);
    %Theta(:,row) = [];
    Xi = Theta\delta;  % initial guess: Least-squares
    
end
