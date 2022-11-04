%% Sparse Identification of Nonlinear Dynamics for discrepancy modeling of the residual error
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
%addpath('./brunton_DMDbook/CODE/CH01_INTRO'); 
addpath('./sparsedynamics - Brunton/utils'); 

%% Generage true dynamics (z) and measurement dynamics (y)

%{
noise = 0.00;
dt = 0.01;
tlength = 50; 
tspan = [0:dt:tlength];
cutoff = 4;
polyorder = 3;
usesine = 0;
lambda = 0.0085;
PD = 0.60 ;  

g = '0.01*y(1).*y(1).*y(1)'; % epsilon discrepancy
% g = ['0.001*y(',num2str(inState),').^3']; % epsilon discrepancy
% g = '0.0';

lowpass_filter = 0; % 1 == filter, 0 == no filter

system = 'Vanderpol';
system = 'Lorenz';
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

%% pool Data  (i.e., build library of nonlinear time series)

Theta = poolData(trainSet.X,n,polyorder,usesine);
m = size(Theta,2);

%% Compute sparse regression for discrepancy g(x): sequential least squares

Xi = sparsifyDynamics(Theta,trainSet.E,lambda,n)

% Xi_lasso = zeros(size(Theta,2),n);
% for i = 1
%     Xi_lasso(:,i) = lasso(Theta,delta,'Lambda',0.2);
% end
%Xi = sparsifyDynamics_simplified(Theta,delta,lambda,n,1)
%Xi = sparsifyDynamics_SSR(Theta,delta)

options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
opts = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1));

%% Reconstruction of Training Data

tspan = trainSet.T;
    
[tR,xR_tmp]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan, zeros(n,1) ,opts);  % approximate

xR = xR_tmp + trainSet.X;

%% Plotting Reconstruction of Training Data

figure, plot(tspan,trainSet.Y,'Linewidth',[2]), hold on, plot(tR,xR,'k--','Linewidth',[2])
legend('True','Learned')
title('Reconstruction of Training Data')
 
%% Forecasting on Test Data

tspan = testSet.T;%(1:500);

switch system
       
    case 'Vanderpol'
        
        %poolDataLIST({'x'},{'y'},Xi,n,polyorder,usesine);

        [tA,xA]=ode45(@(t,z)vanderpol_discrep(t,z,mu,g),testSet.T,trainSet.Y(end,:),options); % true model in test domain
        [tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),testSet.T,zeros(n,1) ,opts);  % approximate
        [tC,xC]=ode45(@(t,x)vanderpol(t,x,mu),testSet.T,trainSet.Y(end,:),options); % ideal model (inadequate)
        
        xA = xA + noise*randn(size(xA));
        xB = xB+xC;
        zoomDisp = 3;
        
    case 'Lorenz'
        
        %poolDataLIST({'x'},{'y'},{'z'},Xi,n,polyorder,usesine);
        
        [tA,xA]=ode45(@(t,x)lorenz_discrep(t,x,sigma,beta,rho,g),tspan,trainSet.Y(end,:),options);   % true model
        [tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,zeros(n,1),opts);  % approximate
        [tC,xC]=ode45(@(t,x)lorenz(t,x,sigma,beta,rho),tspan,trainSet.Y(end,:),options); % ideal model (inadequate)

        xA = xA + noise*randn(size(xA));
        xB = xB+xC;
        zoomDisp = 10;
        
       
end

%% Plotting forecast on test data

figure,
subplot(1,3,1)
plot(tC,xA(:,1)-xC(:,1),'k','Linewidth',[2]), hold on
plot(tC,xA(:,1)-xB(:,1),'r--','Linewidth',[2]), hold on
grid on,
title('Forecasted Discrepancy')
legend('Actual','SINDy')
hold on,

subplot(1,3,2)
plot(tA,xA(:,1),'k','Linewidth',[2]), hold on % Truth
plot(tC,xC(:,1),'b--','Linewidth',[2]), hold on % Plato
plot(tC,xB(:,1),'r--','Linewidth',[2]), hold on % Augmented
grid on,
title('Forecasted Trajectory using Error Discrepancy Model')
legend('True','Plato','Augmented')

subplot(1,3,3)
plot(tA,xA(:,1),'k','Linewidth',[2]), hold on % Truth
plot(tC,xC(:,1),'b--','Linewidth',[2]), hold on % Plato
plot(tC,xB(:,1),'r--','Linewidth',[2]), hold on % Augmented
xlim([tC(1),tC(ceil(length(tC)/zoomDisp))])
grid on,
title('Zoomed in Forecasted Trajectory using Discrepancy Error Model')
legend('True','Plato','Augmented')

sgtitle('Augmenting Platonic Model with Corrective Discrepancy Model of Systematic Error')
hold off

set(gcf,'position',[100,300,1200,400])

%% Calculate RSME of state space error

% clear RMSE_plato RMSE_aug perChange
% 
% RMSE_plato = sqrt(mean((xA(:,1) - xC(:,1)).^2));
% RMSE_aug = sqrt(mean((xA(:,1) - ( xB+xC(:,1) )).^2));
% perChange = (RMSE_plato-RMSE_aug)/(RMSE_plato)*100;
% 
% disp(['RMSE of Platonic state space error: ',num2str(RMSE_plato)])
% disp(['RMSE of Augemented state space error: ',num2str(RMSE_aug)])
% disp(['Percent decrease in RMSE: ',num2str(perChange),'%'])

return;