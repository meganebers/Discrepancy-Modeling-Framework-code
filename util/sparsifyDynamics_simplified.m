function [Xi] = sparsifyDynamics_simplified(Theta,dXdt,lambda,N,NormalizeLib)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data:
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

% compute Sparse regression: sequential least squares
% SINDy with library normalization
% Modified By: K.Kahirman
% Last Updated: 2019/05/14

%% Normalize the library data
if NormalizeLib==1
    for norm_k=1:size(Theta,2)
        normLib(norm_k) = norm(Theta(:,norm_k));
        Theta(:,norm_k) = Theta(:,norm_k)/normLib(norm_k);
    end
end

%% Peform sparse regression
Xi = Theta\dXdt;  % initial guess: Least-squares
[n,m]=size(dXdt);

% lambda is our sparsification knob.
for k=1:N
    smallinds = (abs(Xi)<lambda);   % find small coefficients
    Xi(smallinds)=0;                % and threshold
    for ind = 1:m                   % n is state dimension
        biginds = ~smallinds(:,ind);
        % Regress dynamics onto remaining terms to find sparse Xi
        Xi(biginds,ind) = Theta(:,biginds)\dXdt(:,ind); 
    end
end

%% Now output the SINDy Identified ODEs

% Now retrive the parameters
if NormalizeLib==1
    for norm_k=1:length(Xi)
        Xi(norm_k,:) = Xi(norm_k,:)/normLib(norm_k);
    end
end

