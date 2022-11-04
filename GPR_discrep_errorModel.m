%% Gaussian Process Regression for discrepancy modeling of residual error

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
% needed to "artificially inflate" state space

%{
% p = 50;
% nOrder = size(trainSet.E, 2); % mm1 = m - 1 time_dynamics = zeros(r, mm1);
% q = length(trainSet.E)-p; % length of each Hankel matrix row
% 
% Y_H=[];
% E_H=[];
% for j=1:p
%     Y_H=[Y_H; trainSet.Y(j:q+j,:)];
%     E_H=[E_H; trainSet.E(j:q+j,1)];
% end
% Y_H; % input; state space measurements
% E_H; % output; discrepancy dynamics
%}

%% Perform Gaussian Process regression on discrepancy signal

clear DMdl

for i = 1:n
    T = [trainSet.X trainSet.E(:,i)];
    tbl = array2table(T);
    
    switch system
        
        case 'Lorenz'
            
            tbl.Properties.VariableNames = {'X1','X2','X3','E'};
            
        case 'Vanderpol'
            
            tbl.Properties.VariableNames = {'X1','X2','E'};
            
    end

    DMdl{i} = fitrgp(tbl,'E','KernelFunction','squaredexponential');%,'ActiveSetSize',100,'FitMethod','sr','PredictMethod','fic');

end

%% Evaluate GPR discrepancy model

figure,
for i = 1:n
    xR_tmp(:,i) = resubPredict(DMdl{i});
    subplot(1,n,i), plot(xR_tmp(:,i),'Linewidth',[2]), hold on, plot(trainSet.E(:,i),'--','Linewidth',[2])
    legend('GP','Train')
    hold on,
end
hold off
    
xR = xR_tmp + trainSet.X;

%% Plotting Reconstruction of Training Data

figure, plot(trainSet.T,trainSet.Y,'Linewidth',[2]), hold on, plot(trainSet.T,xR,'k--','Linewidth',[2])
legend('True','Learned')
title('Reconstruction of Training Data')

%% combine GP model and ODE model for forecasting

tspan = testSet.T;
x0 = trainSet.Y(end,:);

switch system
    
   case 'Lorenz'
        
        [tA,xA]=ode45(@(t,x)lorenz_discrep(t,x,sigma,beta,rho,g),tspan,x0,options);   % true model
        xA = xA + noise*randn(size(xA));
       
        %for i = 1:n
        %    xA_butter(:,i) = butterfilt( xA(:,i),cutoff,4, 'Low',ty,'on'  );
        %end
        
        [tC,xC]=ode45(@(t,x)lorenz(t,x,sigma,beta,rho),tspan,x0,options); % ideal model (inadequate)
        
        xB_tmp = [];
        for i = 1:n
            xB_tmp(:,i) = predict(DMdl{i},xC);
        end
        xB = xB_tmp + xC; 

        zoomDisp = 5;
        
    case 'Vanderpol'
        
        [tA,xA]=ode45(@(t,x)vanderpol_discrep(t,x,mu,g),tspan,x0,options);   % true model
        xA = xA + noise*randn(size(xA)); 

        [tC,xC]=ode45(@(t,x)vanderpol(t,x,mu),tspan,x0,options); % ideal model (inadequate)

        xB_tmp = [];
        for i = 1:n
            xB_tmp(:,i) = predict(DMdl{i},xC);
        end
        xB = xB_tmp + xC;
        
        zoomDisp = 3;

end

%% Show augmented model captures true dynamics

figure,
subplot(1,3,1)
plot(tC,xA(:,1) - xC(:,1),'k','Linewidth',[2]), hold on
plot(tC,xA(:,1) - xB(:,1),'r--','Linewidth',[2]), hold on
grid on,
title('Forecasted Discrepancy')
legend('Actual','DMD')
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
% RMSE_aug = sqrt(mean((xA(:,1) - ( pred_approx+xC(:,1) )).^2));
% perChange = (RMSE_plato-RMSE_aug)/(RMSE_plato)*100;
% 
% disp(['RMSE of Platonic state space error: ',num2str(RMSE_plato)])
% disp(['RMSE of Augmented state space error: ',num2str(RMSE_aug)])
% disp(['Percent decrease in RMSE: ',num2str(perChange),'%'])

return;

