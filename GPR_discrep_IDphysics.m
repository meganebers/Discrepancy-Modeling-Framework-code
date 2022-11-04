% Gaussian Process Regression for discrepancy modeling

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
cutoff = 6;
order = 4;
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

%% Evaluate GPR discrepancy model

figure,
for i = 1:n
    xR_tmp(:,i) = resubPredict(DMdl{i});
    subplot(1,n,i), plot(trainSet.Ef(:,i),'Linewidth',[2]), hold on, plot(xR_tmp(:,i),'k--','Linewidth',[2]), hold on,
end
legend('True','Learned')
hold off

%% Reconstruct Training Data

tspan = trainSet.T;
options = odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,n));

switch system
    
    case 'Lorenz'
        
        [tR,xR]=ode45(@(t,x) lorenz_approx_GPR(t,x,sigma,beta,rho,DMdl),tspan,x0,options);  
        
    case 'Vanderpol'
        
        [tR,xR]=ode45(@(t,x) vanderpol_approx_GPR_discrep(t,x,mu,DMdl),tspan,x0,options);  
        
end

%% Plotting Reconstruction of Training Data

figure, plot(trainSet.T,trainSet.Y,'Linewidth',[2]), hold on, plot(trainSet.T,xR,'k--','Linewidth',[2])
legend('True','Learned')
title('Reconstruction of Training Data')

%% combine GP model and ODE model for forecasting

tspan = testSet.T;
x0 = trainSet.Y(end,:);
options = odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,n));

switch system
    
    case 'Lorenz'
        
        [tA,xA]=ode45(@(t,x)lorenz_discrep(t,x,sigma,beta,rho,g),tspan,x0,options);   % true model
        xA = xA + noise*randn(size(xA));
        
        %for i = 1:n
        %    xA_butter(:,i) = butterfilt( xA(:,i),cutoff,4, 'Low',ty,'on'  );
        %end
        
        [tB,xB]=ode45(@(t,x) lorenz_approx_GPR(t,x,sigma,beta,rho,DMdl),tspan,x0,options);  % ideal + discrep model in test domain
        [tC,xC]=ode45(@(t,x) lorenz(t,x,sigma,beta,rho),tspan,x0,options); % ideal model (inadequate)
        
    case 'Vanderpol'
        
        [tA,xA]=ode45(@(t,z) vanderpol_discrep(t,z,mu,g),tspan,x0,options); % true model in test domain
        xA = xA + noise*randn(size(xA));
        
        %for i = 1:n
        %    xA_butter(:,i) = butterfilt( xA(:,i),cutoff,4, 'Low',ty,'on'  );
        %end
        
        [tB,xB]=ode45(@(t,x) vanderpol_approx_GPR_discrep(t,x,mu,DMdl),tspan,x0,options);  % ideal + discrep model in test domain
        [tC,xC]=ode45(@(t,x) vanderpol(t,x,mu),tspan,x0,options); % ideal model (inadequate)
        
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

%% Calculate RSME of dynamical error

% clear RMSE_plato RMSE_aug perChange
% 
% for i = 1:n
%     RMSE_plato(:,i) = sqrt(mean((xA(:,i) - xC(:,i)).^2));
%     RMSE_aug(:,i) = sqrt(mean((xA(:,i) - xB(:,i)).^2));
%     perChange(i) = (RMSE_plato(:,i)-RMSE_aug(:,i))/(RMSE_plato(:,i))*100;
% end
% 
% disp(['RMSE of Platonic dynamical error: ',num2str(RMSE_plato)])
% disp(['RMSE of Augemented dynamical error: ',num2str(RMSE_aug)])
% disp(['Percent decrease in RMSE: ',num2str(perChange),'%'])

return;