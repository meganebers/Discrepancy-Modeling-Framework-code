%% Neural Network for discrepancy modeling of residual error

% Generate data from ideal toy model
% Generate toy truth model (ideal + eps*discrepancy) data
% Generate toy measurement (ideal + eps*discrepancy + noise) data
% Partition data into training and test sets
% Clean noisy signal, if applicable
% Subtract signals to isolate discrepancy
% Perform regression on training discrepancy signal
% Augment ideal model with identified discrepancy regression model
% Forecast augmented ideal model and evaluate error
% Perform data assimilation using new measurement data

%clear all, close all, clc
addpath('./DYNAMICS'); % system dynamics
addpath('./util'); % other functions
addpath('./ODE_Solvers'); % system dynamics
addpath('./sparsedynamics - Brunton/utils');

%% Generage true dynamics (z) and measurement dynamics (y)

%{
noise = 0.000;
dt = 0.01;
tlength = 50;
tspan = [0:dt:tlength];
cutoff = 4;
order = 4;
PD = 0.60 ;  

g = '0.01*y(1).*y(1).*y(1)'; % epsilon discrepancy
% g = ['0.001*y(',num2str(inState),').^3']; % epsilon discrepancy
% g = '0.0';

lowpass_filter = 0; % 1 == filter, 0 == no filter

system = 'Vanderpol'; r = 15; p = 100; % = rank, p = time delays
%system = 'Vanderpol'; r = 50; % = rank, p = time delays
%}

eval(['discrepancyDynamics_',system])

T = tx;
E = e; % state space error
X = x;
Y = y;
Z = z;

trainSet.T = T(1:round(PD*length(T)),:);
testSet.T = T(round(PD*length(T)):end,:);

trainSet.E = E(1:round(PD*length(E)),:);
testSet.E = E(round(PD*length(E)):end,:);

trainSet.X = X(1:round(PD*length(X)),:);
testSet.X = X(round(PD*length(X)):end,:);

trainSet.Y = Y(1:round(PD*length(Y)),:);
testSet.Y = Y(round(PD*length(Y)):end,:);

trainSet.Z = Z(1:round(PD*length(Z)),:);
testSet.Z = Z(round(PD*length(Z)):end,:);

%% Time Delay for continuous time system
% to "artificially inflate" state space

%{
p = 50;
nOrder = size(trainSet.E, 2); % mm1 = m - 1 time_dynamics = zeros(r, mm1);
q = length(trainSet.E)-p; % length of each Hankel matrix row
 
Y_H=[];
E_H=[];
for j=1:p
    Y_H=[Y_H; trainSet.X(j:q+j,:)];
    E_H=[E_H; trainSet.E(j:q+j,1)];
end
Y_H; % input; state space measurements
E_H; % output; discrepancy dynamics
%}

%% Test NN residual model

%
clear net
numNodes = [10 10 10];
net = feedforwardnet(numNodes,'trainbr');
net.trainParam.epochs = 1000;
net.layers{1}.transferFcn = 'logsig';
net.layers{2}.transferFcn = 'radbas';
net.layers{3}.transferFcn = 'purelin';
%net = train(net,Y_H.',E_H.','useParallel','yes','useGPU','only');
net = train(net,trainSet.X.',trainSet.E.')%,'useParallel','yes','useGPU','only');
yy = net(trainSet.X.');
perf_train = perform(net,yy,trainSet.E.')

figure, 
subplot(1,2,1), plot(yy.'), hold on, plot(trainSet.E,'--')

EPred=[];
for i = 1:length(testSet.E)
    [EPred(i,:)] = net(testSet.X(i,:).');
end

subplot(1,2,2), plot(EPred), hold on, plot(testSet.E,'--')

perf_test = perform(net,EPred.',testSet.E.')

xR = yy.' + trainSet.X;

%% Plotting Reconstruction of Training Data

figure, plot(trainSet.T,trainSet.Y,'Linewidth',[2]), hold on, plot(trainSet.T,xR,'k--','Linewidth',[2])
legend('True','Learned')
title('Reconstruction of Training Data')

%% combine GP model and ODE model to check if discrepancy is captured (test domain)

tspan = testSet.T;
x0 = trainSet.Y(end,:);

switch system
    
    case 'Lorenz'
        
        [tA,xA]=ode45(@(t,x)lorenz_discrep(t,x,sigma,beta,rho,g),tspan,x0,options);   % true model in test domain
        xA = xA + noise*randn(size(xA));
        
%         xA_butter = [];
%         %cutoff = ceil(std(xA));
%         for i = 1:n
%             xA_butter(:,i) = butterfilt( xA(:,i),cutoff,order, 'Low',ty,'on'  );
%         end
        
        [tC,xC]=ode45(@(t,x)lorenz(t,x,sigma,beta,rho),tspan,x0,options); % ideal model (inadequate)
       
        xB_tmp=[];
        for i = 1:length(xC)
            [xB_tmp(i,:)] = net(xC(i,:).');
        end
        xB = xB_tmp + xC;

    case 'Vanderpol'
        
        [tA,xA]=ode45(@(t,x)vanderpol_discrep(t,x,mu,g),tspan,x0,options);   % true model in test domain
        xA = xA + noise*randn(size(xA));
        
%         xA_butter = [];
%         %cutoff = ceil(std(xA));
%         for i = 1:n
%             xA_butter(:,i) = butterfilt( xA(:,i),cutoff,order, 'Low',ty,'on'  );
%         end
        
        [tC,xC]=ode45(@(t,x)vanderpol(t,x,mu),tspan,x0,options); % ideal model (inadequate)
        
        xB_tmp=[];
        for i = 1:length(xC)
            [xB_tmp(i,:)] = net(xC(i,:).');
        end
        xB = xB_tmp + xC;
end

%% Show augmented model captures true dynamics

figure,
subplot(1,2,1), plot(tC, xA(:,1)-xC(:,1),'Linewidth',[2]), hold on, plot(tC,xA(:,1)-xB(:,1),'--','Linewidth',[2])
title('Discrepancy dynamics')

subplot(1,2,2)
plot(tA,xA(:,1),'k','Linewidth',[2]), hold on
plot(tC,xB(:,1),'r--','Linewidth',[2]), hold on
plot(tC,xC(:,1),'b--','Linewidth',[2]), hold on
grid on,
title('Forecasted Trajectory using Discrepancy Model')
legend('True','Augmented','Plato')

return;