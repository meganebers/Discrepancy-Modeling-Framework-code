function [Phi,omega,lambda,b,Xdmd,D,Atilde,U_r,V_r,W_r] = DMD_simplified(X1,X2,r,dt)
% function [Phi,omega,lambda,b,Xdmd] = DMD(X1,X2,r,dt)
% Computes the Dynamic Mode Decomposition of X1, X2
%
% INPUTS:
% X1 = X, data matrix
% X2 = X', shifted data matrix
% Columns of X1 and X2 are state snapshots
% r = target rank of SVD
% dt = time step advancing X1 to X2 (X to X')
%
% OUTPUTS:
% Phi, the DMD modes
% omega, the continuous-time DMD eigenvalues
% lambda, the discrete-time DMD eigenvalues
% b, a vector of magnitudes of modes Phi
% Xdmd, the data matrix reconstrcted by Phi, omega, b

%% DMD
[U, S, V] = svd(X1, 'econ');
r = min(r, size(U,2));
if r == size(U,2)
    a = 1;
else
    a = 2; % set scale for visualizing modes
end

figure,
subplot(2,1,1), plot(diag(S)/(sum(diag(S))),'ro','Linewidth',[3]), title('SVD Modes')
subplot(2,1,2), plot(diag(S(1:r*a,1:r*a))/(sum(diag(S(1:r*a,1:r*a)))),'ro','Linewidth',[3]), title('Truncated SVD Modes')

figure, subplot(2,1,1), plot(U(:,1:r),'Linewidth',[2]), title('SVD Modes')
subplot(2,1,2), plot(V(:,1:r),'Linewidth',[2]), title('SVD Modes')

U_r = U(:, 1:r); % truncate to rank-r
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);
% Atilde = U_r' * X2 * V_r / S_r; % low-rank dynamics


Atilde = U_r' * X2 * V_r / S_r;
[W_r, D] = eig(Atilde);
Phi = X2 * V_r / S_r * W_r; % DMD modes


lambda = diag(D); % discrete-time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues

% Check attractor structure
% figure,
% plot(V_r(:,1),V_r(:,2))
% 
% figure,
% plot3(V_r(:,1),V_r(:,2),V_r(:,3))

%% Compute DMD mode amplitudes b
x1 = X1(:, 1);
b = Phi\x1;

%% DMD reconstruction
mm1 = size(X1, 2); % mm1 = m - 1
time_dynamics = zeros(r, mm1);
t = (0:mm1-1)*dt; % time vector

for iter = 1:mm1
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end

Xdmd = Phi * time_dynamics;
% XHat = Xdmd(1:2,:).';
% 
% figure,
% plot(XHat)
% 
% for j=1:idx
%   H=[H; Dtrain(j:idx+j,:).'];
% end  




