%% Sparse Identification of Nonlinear Dynamics for discrepancy modeling & learning missing physics

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

%clear all, close all, clc;
%addpath('./Method Scripts'); % system dynamics
addpath('./DYNAMICS'); % system dynamics
addpath('./util'); % other functions
%addpath('./DYNAMICS/PDE'); % system dynamics
addpath('./sparsedynamics - Brunton/utils'); 

%% Generage true dynamics (z) and measurement dynamics (y)

%
noise = 0.00;
dt = 0.001;
tlength = 50; 
tspan = [0:dt:tlength];
cutoff = 4;
polyorder = 3;
usesine = 0;
lambda = 0.0095; % lambda is our sparsification knob; Vanderpol (0.0085)
PD = 0.60 ;  

g = '0.01*y(1).*y(1).*y(1)'; % epsilon discrepancy
% g = ['0.001*y(',num2str(inState),').^3']; % epsilon discrepancy
% g = '0.0';

lowpass_filter = 0; % 1 == filter, 0 == no filter

system = 'Vanderpol';
%}

eval(['discrepancyDynamics_',system])

T = tx;
Ef = ef; % dynamical model error
X = x;
Y = y;
Z = z;

trainSet.T = T(1:round(PD*length(T))); 
testSet.T = T(round(PD*length(T)):end);

trainSet.Ef = Ef(1:round(PD*length(Ef)),:);
testSet.Ef = Ef(round(PD*length(Ef)):end,:);

trainSet.X = X(1:round(PD*length(X)),:);
testSet.X = X(round(PD*length(X)):end,:);

trainSet.Y = Y(1:round(PD*length(Y)),:);
testSet.Y = Y(round(PD*length(Y)):end,:);

trainSet.Z = Z(1:round(PD*length(Z)),:);
testSet.Z = Z(round(PD*length(Z)):end,:);

%% pool Data  (i.e., build library of nonlinear time series)

Theta = poolData(trainSet.Y,n,polyorder,usesine);
m = size(Theta,2);

%% Compute sparse regression for discrepancy g(x): sequential least squares

Xi = sparsifyDynamics(Theta,trainSet.Ef,lambda,n)

% Xi_lasso = zeros(size(Theta,2),n);
% for i = 1
%     Xi_lasso(:,i) = lasso(Theta,delta,'Lambda',0.2);
% end
%Xi = sparsifyDynamics_simplified(Theta,delta,lambda,n,1)
%Xi = sparsifyDynamics_SSR(Theta,trainSet.Ef)

options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));

%% Reconstruction of Training Data

tspan = trainSet.T;
    
switch system
    
    case 'Vanderpol'
        
        [tR,xR]=ode45(@(t,x)vanderpol_approx_SINDy_discrep(t,x,Xi,mu,polyorder,usesine),tspan,x0,options);
        
    case 'Lorenz'
        
        [tR,xR]=ode45(@(t,x)lorenz_approx_SINDy_discrep(t,x,Xi,sigma,beta,rho,polyorder,usesine),tspan,x0,options);
        
  
end

%% Plotting Reconstruction of Training Data

figure, plot(tspan,trainSet.Y,'Linewidth',[2]), hold on, plot(tR,xR,'k--','Linewidth',[2])
legend('True','Learned')
title('Reconstruction of Training Data')

%% Forecasting on Test Data
switch system
         
    case 'Vanderpol'
        
        poolDataLIST({'x','y'},Xi,n,polyorder,usesine);
        
        [tA,xA]=ode45(@(t,z)vanderpol_discrep(t,z,mu,g),testSet.T,trainSet.Y(end,:),options); % true model in test domain
        [tB,xB]=ode45(@(t,x)vanderpol_approx_SINDy_discrep(t,x,Xi,mu,polyorder,usesine),testSet.T,trainSet.Y(end,:),options);  % ideal + discrep model in test domain
        [tC,xC]=ode45(@(t,x)vanderpol(t,x,mu),testSet.T,trainSet.Y(end,:),options); % ideal model (inadequate)
           
        xA = xA + noise*randn(size(xA));
        
    case 'Lorenz'
        
        poolDataLIST({'x','y','z'},Xi,n,polyorder,usesine);
        
        [tA,xA]=ode45(@(t,z)lorenz_discrep(t,z,sigma,beta,rho,g),testSet.T,trainSet.Y(end,:),options);   % true model
        [tB,xB]=ode45(@(t,x)lorenz_approx_SINDy_discrep(t,x,Xi,sigma,beta,rho,polyorder,usesine),testSet.T,trainSet.Y(end,:),options);  % approximate
        [tC,xC]=ode45(@(t,x)lorenz(t,x,sigma,beta,rho),testSet.T,trainSet.Y(end,:),options); % ideal model (inadequate)

        xA = xA + noise*randn(size(xA));
        
       
end

%% Plotting Forecasting Test Data

figure,
for i = 1:n
    subplot(n,1,i)
    plot(tA,xA(:,i),'k','Linewidth',[2]), hold on % Truth
    plot(tC,xC(:,i),'b--','Linewidth',[2]), hold on % Plato
    plot(tB,xB(:,i),'r--','Linewidth',[2]), hold on % Augmented
    xlim([tC(1),tC(roundn(length(tC), 2)/2)])
    grid on,
end
sgtitle('Forecasting Over Training Data')
%legend('True','Plato','Augmented') 

set(gcf,'position',[100,300,1200,400],'DefaultFigureRenderer', 'painters')

%%

state = 1;

figure, 
plot(tspan,trainSet.X(:,state),'b-','Linewidth',[5])
hold on
plot(tspan,trainSet.Y(:,state),'r-','Linewidth',[5]), 
hold on,
plot(tR,xR(:,state),'k--','Linewidth',[5]),
hold on, 
plot(tC,xC(:,state),'b-','Linewidth',[5]), hold on % Plato
hold on
plot(tB,xB(:,state),'r-','Linewidth',[5]),
hold on, 
plot(tA,xA(:,state),'k--','Linewidth',[5]), 
hold on,
xline(tA(1),'k','Linewidth',[2])

state = 2;

figure, 
plot(tspan,trainSet.X(:,state),'b-','Linewidth',[5])
hold on
plot(tspan,trainSet.Y(:,state),'r-','Linewidth',[5]), 
hold on,
plot(tR,xR(:,state),'k--','Linewidth',[5]),
hold on, 
plot(tC,xC(:,state),'b-','Linewidth',[5]), hold on % Plato
hold on
plot(tB,xB(:,state),'r-','Linewidth',[5]),
hold on, 
plot(tA,xA(:,state),'k--','Linewidth',[5]), 
hold on,
xline(tA(1),'k','Linewidth',[2])


return;
