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
%addpath('./brunton_DMDbook/CODE/CH01_INTEO'); 
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
%system = 'Vanderpol'; r = 50; p = 300; % = rank, p = time delays
 %}

eval(['discrepancyDynamics_',system])

T = tx;
Ef = ef; % state space error
X = x;
Y = y;
Z = z;

trainSet.T = T(1:round(PD*length(T)),:); 
testSet.T = T(round(PD*length(T)):end,:);

trainSet.Ef = Ef(1:round(PD*length(Ef)),:);
testSet.Ef = Ef(round(PD*length(Ef)):end,:);

trainSet.X = X(1:round(PD*length(X)),:);
testSet.X = X(round(PD*length(X)):end,:);

trainSet.Y = Y(1:round(PD*length(Y)),:);
testSet.Y = Y(round(PD*length(Y)):end,:);

trainSet.Z = Z(1:round(PD*length(Z)),:);
testSet.Z = Z(round(PD*length(Z)):end,:);

%% Train Neural Network

clear net
numNodes = [10 10 10];
net = feedforwardnet(numNodes,'trainbr');
net.trainParam.epochs = 10000;
net.layers{1}.transferFcn = 'logsig';
net.layers{2}.transferFcn = 'radbas';
net.layers{3}.transferFcn = 'purelin';
%net = fitnet(numNodes,'trainbr');
net = train(net,trainSet.Y.',trainSet.Ef.');%,'useParallel','yes','useGPU','only');
yy = net(trainSet.Y.');
perf = perform(net,yy,trainSet.Ef.')

%% Reconstruct Training Data

tspan = trainSet.T;
        
switch system
    
    case 'Lorenz'
        
        [tR,xR]=ode45(@(t,x)lorenz_approx_NN_discrep(t,x,sigma,beta,rho,net),tspan,x0,options); % augmented model
                
    case 'Vanderpol'
        
        [tR,xR]=ode45(@(t,x)vanderpol_approx_NN_discrep(t,x,mu,net),tspan,x0,options); % augmented model
            
end

%% Plotting Reconstruction of Training Data

figure, plot(trainSet.T,trainSet.Y,'Linewidth',[2]), hold on, plot(trainSet.T,xR,'k--','Linewidth',[2])
legend('True','Learned')
title('Reconstruction of Training Data')

%% combine NN model and ODE model to check if discrepancy is captured (test domain)

x0 = trainSet.Y(end,:);
tspan = testSet.T;
        
switch system
    
    case 'Lorenz'
        
        [tA,xA]=ode45(@(t,x)lorenz_discrep(t,x,sigma,beta,rho,g),tspan,x0,options);   % true model in test domain
        xA = xA + noise*randn(size(xA));
        
        [tB,xB]=ode45(@(t,x)lorenz_approx_NN_discrep(t,x,sigma,beta,rho,net),tspan,x0,options); % augmented model

        [tC,xC]=ode45(@(t,x)lorenz(t,x,sigma,beta,rho),tspan,x0,options); % ideal model (inadequate)
                
    case 'Vanderpol'
        
        [tA,xA]=ode45(@(t,x)vanderpol_discrep(t,x,mu,g),tspan,x0,options);   % true model in test domain
        xA = xA + noise*randn(size(xA));
        
        [tB,xB]=ode45(@(t,x)vanderpol_approx_NN_discrep(t,x,mu,net),tspan,x0,options); % augmented model

        [tC,xC]=ode45(@(t,x)vanderpol(t,x,mu),tspan,x0,options); % ideal model (inadequate)
            
end

%% Plotting Results

figure,
for i = 1:n
    subplot(1,n,i)
    plot(tA,xA(:,i),'k','Linewidth',[2]), hold on % Truth
    plot(tC,xC(:,i),'b--','Linewidth',[2]), hold on % Plato
    plot(tB,xB(:,i),'r--','Linewidth',[2]), hold on % Augmented
    grid on,
end
sgtitle('Forecasted Trajectory using Discrepancy Model')
legend('True','Plato','Augmented')

figure,
for i = 1:n
    subplot(1,n,i)
    plot(tA,xA(:,i),'k','Linewidth',[2]), hold on % Truth
    plot(tC,xC(:,i),'b--','Linewidth',[2]), hold on % Plato
    plot(tB,xB(:,i),'r--','Linewidth',[2]), hold on % Augmented
    xlim([tC(1),tC(roundn(length(tC), 2)/2)])
    grid on,
end
sgtitle('Forecasted Trajectory using Discrepancy Model')
legend('True','Plato','Augmented')

set(gcf,'position',[100,300,1200,400])

return;