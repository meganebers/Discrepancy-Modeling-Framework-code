%%
clear all; close all; clc

addpath('./sparsedynamics - Brunton/utils')
addpath('./optdmd-master/src'); % optDMD
addpath('./DYNAMICS'); 

%% Get state trajectories
 
noise = 0.0;
eval(['Mass_Spring_Damper'])

y_M = y_L(:,1);
y_E = y_L_bias; % Apply measurement bias to the output

E = y_M(:,1)-y_E(:,1);

% Plot results
figure, subplot(1,2,1),
plot(t, y_M(:,1), '-', t, y_E, '--', 'LineWidth', 2);
legend({'Position', 'Position with Measurement Bias'});
ylabel('Position [m]')
xlabel('Time [s]')
title('Mass-Spring-Damper system with Systematic Residual');
grid on;
hold on

% Plot results
subplot(1,2,2),
plot(t, E, '-', 'LineWidth', 2);
legend({'Position Error'});
ylabel('Position Error')
xlabel('Time [s]')
title('Error of Mass-Spring-Damper system vs Measurement');
grid on;
hold off

%%
%{
% System parameters
m = 2;  % Mass
k = 2;  % Spring constant
b = 1;  % Damping coefficient

% State-space matrices
A = [0 1; -k/m -b/m];
B = [0; 1/m];
C = [1 0];
C_bias = [1+bias 0];
D = 0;

% Create state-space model
sys = ss(A, B, C, D);
sys_bias = ss(A, B, C_bias, D);

% Simulate the system
[y_M, ~, x] = lsim(sys, zeros(size(t)), t, S);
[y_E, ~, x_bias] = lsim(sys_bias, zeros(size(t)), t, S);

plot(t, y_M, '-', t, y_E, '--', 'LineWidth', 2);
legend({'Position', 'Position with Measurement Bias'});
ylabel('Position [m]')
xlabel('Time [s]')
title('Mass-Spring-Damper system with Measurement Bias');
hold off
%}

%% SINDy
m = size(E,2);
Theta = poolData(y_M,m,5,0);

lambda = 0.025; % lambda is our sparsification knob; Vanderpol (0.0085)
Xi = sparsifyDynamics(Theta,E,lambda,m)

y_tilde_SINDy = y_E+Theta*Xi;

figure(2), subplot(2,2,1)
plot(t, y_tilde_SINDy, 'k-', t, y_M, 'r--', 'Linewidth',[3])
hold on, title('SINDy ')
xlim([0, 20]),
grid on

%% optDMD

% Time Delay for continuous time system
% needed to "artificially inflate" state space
p = 10;
r = 5;
nOrder = size(E, 2); % mm1 = m - 1 time_dynamics = zeros(r, mm1);
q = length(E)-p; % length of each Hankel matrix row

Y_H=[];
D_H=[];
for j=1:p
    %Y_H=[Y_H; trainSet.Y(j:q+j,:).'];
    D_H=[D_H; E(j:q+j,:).'];
end
Y_H; % input; state space measurements
D_H; % output; discrepancy dynamics

% optDMD
% linear constraints; meant to constrain eigs to neg real
lbc = [-Inf*ones(r); -Inf*ones(r)];
ubc = [zeros(r); Inf*ones(r)];

copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);

ts = t(1:q+1);
dt = ts(2)-ts(1);
imode = 3;
[DMD.Phi,DMD.lambda,DMD.b] = optdmd(D_H,ts,r,imode,[],[],[],copts);
DMD.omega = log(DMD.lambda)/dt; % discrete to continuous time eigs
disp('DMD computed')

% Reconstruction

%D_recon = w*diag(b)*exp(e*t')
D1 = DMD.Phi*diag(DMD.b)*exp(DMD.lambda*t.');
%relerr_r = norm(D1-Y_H,'fro')/norm(Y_H,'fro')
%relerr_r_clean = norm(D1-xclean,'fro')/norm(xclean,'fro');

xR_tmp = real(D1(1:nOrder,:)).';

figure,
plot(D_H(1:nOrder,:).','Linewidth',[2]), hold on, % cols is
plot(xR_tmp,'--','Linewidth',[2]), hold on,
plot(E)
legend('True','','','Reconstructed','','')
title('DMD Reconstruction of First Time Delay')

y_tilde_DMD = y_E + xR_tmp;

figure(2), subplot(2,2,2)
plot(t, y_tilde_DMD, 'k-', t, y_M, 'r--', 'Linewidth',[3])
hold on, title('optDMD')
xlim([0, 20]),
grid on

%% GPR

clear DMdl

T = [y_M E];
tbl = array2table(T);

tbl.Properties.VariableNames = {'y_M','E'};

DMdl = fitrgp(tbl,'E','KernelFunction','squaredexponential');%,'ActiveSetSize',100,'FitMethod','sr','PredictMethod','fic');

xR_tmp = resubPredict(DMdl);

y_tilde_GPR = y_E + xR_tmp;

figure(2), subplot(2,2,3)
plot(t, y_tilde_GPR, 'k-', t, y_M, 'r--', 'Linewidth',[3])
xlim([0, 20]),
grid on
hold on, title('GPR')

%% NN
clear net
numNodes = [10 10 10];
net = feedforwardnet(numNodes,'trainbr');
net.trainParam.epochs = 500;
net.layers{1}.transferFcn = 'logsig';
net.layers{2}.transferFcn = 'radbas';
net.layers{3}.transferFcn = 'purelin';
%net = train(net,Y_H.',E_H.','useParallel','yes','useGPU','only');
net = train(net,y_E',E');%,'useParallel','yes','useGPU','only');
%yy = net(y_E');
%perf_train = perform(net,yy,trainSet.E.')

yy_tilde_NN = y_E + net(y_E').';

figure(2), subplot(2,2,4)
plot(t, yy_tilde_NN, 'k-', t, y_M, 'r--', 'Linewidth',[3])
hold on, title('NN')
xlim([0, 20]),
grid on

%

figure(2), sgtitle('Observation Error of Simple Moving Object')
set(gcf, 'PaperType', 'usletter');