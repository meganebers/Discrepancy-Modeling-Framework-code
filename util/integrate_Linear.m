function [dxdt] = integrate_Linear(t, X, U, ut, A, B, Tfinal)
%integrate_Lineear Integrate a linear model using ode integrator
%
% Inputs
% t - time
% X - states [state x time ]
% A - [output x state] state transition matrix
% B - [output x input] input matrix
% U - [ input x time] inputs
% ut - input time vecotr
% Tfinal - final time

% M Ebers, modified from M Rosenberg - UW - 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%% 1. Identify hybrid region

phs = round(t/Tfinal*100);
for ii = 1:length(phs)
%% 2. Compute derivatives
    if isempty(B)
        dxdt(:,ii) = A*X(:,ii);
    else
        u = interp1(ut, U',t)';
        dxdt(:,ii) = A*X(:,ii) + B*u(:,ii);
    end
end

end

