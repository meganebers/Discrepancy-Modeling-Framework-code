%% Dynamic Mode Decomposition for discrepancy modeling of residual error

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
%set(0,'DefaultFigureVisible','on')
addpath('./DYNAMICS'); % system dynamics
addpath('./DYNAMICS/PDE'); % system dynamics
addpath('./util'); % other functions
addpath('./optdmd-master/src'); % optDMD
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

system = 'Vanderpol'; r = 15; p = 100; % = rank, p = time delays
system = 'Lorenz'; r = 50; p = 300; % = rank, p = time delays
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

nOrder = size(trainSet.E, 2); % mm1 = m - 1 time_dynamics = zeros(r, mm1);
q = length(trainSet.E)-p; % length of each Hankel matrix row

Y_H=[];
D_H=[];
for j=1:p
    %Y_H=[Y_H; trainSet.Y(j:q+j,:).'];
    D_H=[D_H; trainSet.E(j:q+j,:).'];
end
Y_H; % input; state space measurements
D_H; % output; discrepancy dynamics

%% Perform optDMD

% linear constraints; meant to constrain eigs to neg real
lbc = [-Inf*ones(r); -Inf*ones(r)];
ubc = [zeros(r); Inf*ones(r)];

copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);

ts = trainSet.T(1:q+1);
imode = 2;
[DMD.Phi,DMD.lambda,DMD.b] = optdmd(D_H,ts,r,imode,[],[],[],copts);
DMD.omega = log(DMD.lambda)/dt; % discrete to continuous time eigs
disp('DMD computed')

% Reconstruction

%D_recon = w*diag(b)*exp(e*t')
D1 = DMD.Phi*diag(DMD.b)*exp(DMD.lambda*trainSet.T.');
%relerr_r = norm(D1-Y_H,'fro')/norm(Y_H,'fro')
%relerr_r_clean = norm(D1-xclean,'fro')/norm(xclean,'fro');

xR_tmp = real(D1(1:nOrder,:)).';

figure,
plot(D_H(1:nOrder,:).','Linewidth',[2]), hold on, % cols is
plot(xR_tmp,'--','Linewidth',[2]), hold off
legend('True','','','Reconstructed','','')
title('DMD Reconstruction of First Time Delay')

xR = xR_tmp + trainSet.X;

%% Plotting Reconstruction of Training Data

figure, plot(trainSet.T,trainSet.Y,'Linewidth',[2]), hold on, plot(trainSet.T,xR,'k--','Linewidth',[2])
legend('True','Learned')
title('Reconstruction of Training Data')

%% combine DMD model and toy system for forecasting

tspan = testSet.T.';
x0 = trainSet.Y(end,:).';
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));

switch system
    
    case 'Lorenz'
        
        [tA,xA]=ode45(@(t,x)lorenz_discrep(t,x,sigma,beta,rho,g),tspan,x0,options);   % true model in test domain
        xA = xA + noise*randn(size(xA));
        
        [tC,xC]=ode45(@(t,x)lorenz(t,x,sigma,beta,rho),tspan,x0,options); % ideal model (inadequate)
        
        zoomDisp = 5;
        
    case 'Vanderpol'
        
        [tA,xA]=ode45(@(t,x)vanderpol_discrep(t,x,mu,g),tspan,x0,options);   % true model in test domain
        xA = xA + noise*randn(size(xA));
        
        [tC,xC]=ode45(@(t,x)vanderpol(t,x,mu),tspan,x0,options); % ideal model (inadequate)
        
        zoomDisp = 3;
        
end

tspan = [0:dt:(length(testSet.T)-1)*dt]; % MUST START AT 0 due to discrepancy error == 0 at initialization regardless of IC, but DMD is dependent on time series and needs to be initialized at t = 0;
xB_tmp = real( DMD.Phi*diag(DMD.b)*exp(DMD.lambda*tspan) );
xB = xB_tmp(1:nOrder,:).' + xC ;

%% Augment ideal model with identified regression discrepancy model

figure,
subplot(1,3,1)
plot(tC,xA(:,1)-xC(:,1),'k','Linewidth',[2]), hold on
plot(tC,xA(:,1)-xB(:,1),'r--','Linewidth',[2]), hold on
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
% RMSE_aug = sqrt(mean((xA(:,1) - ( xD(:,1)+xC(:,1) )).^2));
% perChange = (RMSE_plato-RMSE_aug)/(RMSE_plato)*100;
% 
% disp(['RMSE of Platonic state space error: ',num2str(RMSE_plato)])
% disp(['RMSE of Augemented state space error: ',num2str(RMSE_aug)])
% disp(['Percent decrease in RMSE: ',num2str(perChange),'%'])

return;