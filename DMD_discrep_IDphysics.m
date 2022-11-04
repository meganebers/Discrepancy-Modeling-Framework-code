 %% Dynamic Mode Decomposition for discrepancy modeling

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
addpath('./DYNAMICS/PDE'); % system dynamics
addpath('./util'); % other functions
addpath('./optdmd-master/src'); % optDMD
addpath('./ODE_Solvers');
addpath('./sparsedynamics - Brunton/utils');

%% Generage true dynamics (z) and measurement dynamics (y)

%{
noise = 0.00;
dt = 0.01;
tlength = 50;
tspan = [0:dt:tlength];
cutoff = 10;
order = 4;
PD = 0.60 ;  

g = '0.01*y(1).*y(1).*y(1)'; % epsilon discrepancy
% g = ['0.001*y(',num2str(inState),').^3']; % epsilon discrepancy
% g = '0.0';

lowpass_filter = 1; % 1 == filter, 0 == no filter

system = 'Vanderpol'; r = 15; p = 100; % = rank, p = time delays
system = 'Lorenz'; r = 35; p = 300; % = rank, p = time delays

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

%% Time Delay for continuous time system
% needed to "artificially inflate" state space

nOrder = size(trainSet.Ef, 2); % mm1 = m - 1 time_dynamics = zeros(r, mm1);
q = length(trainSet.Ef)-p; % length of each Hankel matrix row

Y_H=[];
D_H=[];
for j=1:p
    %Y_H=[Y_H; trainSet.Y(j:q+j,:).'];
    D_H=[D_H; trainSet.Ef(j:q+j,:).'];
end
Y_H; % input; state space measurements
D_H; % output; discrepancy dynamics

%% Perform optDMD

% linear constraints; meant to constrain eigs to neg real
lbc = [-Inf*ones(r); -Inf*ones(r)];
ubc = [zeros(r); Inf*ones(r)];

copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);

ts = trainSet.T(1:q+1);
imode = 1;
[DMD.Phi,DMD.lambda,DMD.b] = optdmd(D_H,ts,r,imode,[],[],[],copts);
DMD.omega = log(DMD.lambda)/dt; % discrete to continuous time eigs
disp('DMD computed')

%% Reconstruction of Training Data

D1 = DMD.Phi*diag(DMD.b)*exp(DMD.lambda*ts.');

figure,
for i = 1:nOrder
    subplot(1,nOrder,i)
    plot(D_H(i,:).','Linewidth',[2]), hold on, % cols is
    plot(real(D1(i,:).'),'--','Linewidth',[2]), hold on
    grid on,
end
legend('True','Reconstructed')
sgtitle('DMD Reconstruction of First Time Delay')
hold off

tspan = trainSet.T;

switch system
    
    case 'Lorenz'
        
        x0new = x0;
        xR = [];
        xR = x0new.';
        for i = 1:length(tspan)-1
            tt = [tspan(i) tspan(i+1)];
            xB_tmp = ode5(@(t,x)lorenz_approx_DMD(t,x,sigma,beta,rho,DMD,dt,tt(1)),tt,x0new.');
            x0new = xB_tmp(2,:);
            xR = [xR; x0new];
        end
        tR = tspan;
        
    case 'Vanderpol'
        
        x0new = x0;
        xR = [];
        xR = x0new.';
        for i = 1:length(tspan)-1
            tt = [tspan(i) tspan(i+1)];
            xB_tmp = ode5(@(t,x)vanderpol_approx_DMD(t,x,mu,DMD,dt,tt(1)),tt,x0new.');
            x0new = xB_tmp(2,:);
            xR = [xR; x0new];
        end
        tR = tspan;
        
end

%% Plotting Reconstruction of Training Data

figure, plot(tR,trainSet.Y,'Linewidth',[2]), hold on, plot(tR,xR,'k--','Linewidth',[2])
legend('True','Learned')
title('Reconstruction of Training Data')

%% combine DMD model and toy system for forecasting

x0 = trainSet.Y(end,:).';
tspan = testSet.T;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));

switch system
    
    case 'Lorenz'
        
        [tA,xA]=ode45(@(t,x)lorenz_discrep(t,x,sigma,beta,rho,g),tspan,x0,options);   % true model in test domain
        xA = xA + noise*randn(size(xA));
        
%         for i = 1:n
%             xA_tmp(:,i) = butterfilt( xA(:,i),cutoff,4, 'Low',ty,'on'  );
%         end
%         xA = xA_tmp; 
        
        x0new = x0;
        xB = [];
        xB = x0new.';
        for i = 1:length(tspan)-1
            tt = [tspan(i) tspan(i+1)];
            xB_tmp = ode5(@(t,x)lorenz_approx_DMD(t,x,sigma,beta,rho,DMD,dt,tt(1)),tt,x0new.');
            x0new = xB_tmp(2,:);
            xB = [xB; x0new];
        end
        tB = tspan;
        
        [tC,xC]=ode45(@(t,x)lorenz(t,x,sigma,beta,rho),tspan,x0,options); % ideal model (inadequate)
        
        zoomDisp = 5;
        
    case 'Vanderpol'
        
        
        [tA,xA]=ode45(@(t,x)vanderpol_discrep(t,x,mu,g),tspan,x0,options);   % true model in test domain
        xA = xA + noise*randn(size(xA));
           
%         if lowpass_filter == 1
%             
%             for i = 1:n
%                 xA_tmp(:,i) = butterfilt( xA(:,i),cutoff,4, 'Low',ty,'off'  );
%             end
%             xA = xA_tmp;
%             
%         else
%             
%             xA = xA; % no filtering
%             
%         end
        
        
        x0new = x0;
        xB = [];
        xB = x0new.';
        for i = 1:length(tspan)-1
            tt = [tspan(i) tspan(i+1)];
            xB_tmp = ode5(@(t,x)vanderpol_approx_DMD(t,x,mu,DMD,dt,tt(1)),tt,x0new.');
            x0new = xB_tmp(2,:);
            xB = [xB; x0new];
        end
        tB = tspan;
        
        [tC,xC]=ode45(@(t,x)vanderpol(t,x,mu),tspan,x0,options); % ideal model (inadequate)
        
        zoomDisp = 3;
        
end

%% Augment ideal model with identified regression discrepancy model

figure,
for i = 1:nOrder
    subplot(1,nOrder,i)
    plot(tA,xA(:,i)-xC(:,i),'k','Linewidth',[2]), hold on
    plot(tA,xA(:,i)-xB(:,i),'r--','Linewidth',[2]), hold on
    grid on,
    title(['State ',num2str(i)])
end
sgtitle('Test Domain: Discrepancy Comparison')
legend('Actual','DMD')
hold off
set(gcf,'position',[100,300,800,400])

figure,
for i = 1:nOrder
    subplot(2,nOrder,i)
    plot(tA,xA(:,i),'k','Linewidth',[2]), hold on % Truth
    plot(tB,xB(:,i),'r--','Linewidth',[2]), hold on % Augmented
    plot(tC,xC(:,i),'b--','Linewidth',[2]), hold on % Plato
    grid on,
    title(['State ',num2str(i)])
    
    subplot(2,nOrder,i+nOrder)
    plot(tA,xA(:,i),'k','Linewidth',[2]), hold on % Truth
    plot(tB,xB(:,i),'r--','Linewidth',[2]), hold on % Augmented
    plot(tC,xC(:,i),'b--','Linewidth',[2]), hold on % Plato
    grid on,
    title(['ZOOMED State ',num2str(i)])
    xlim([tC(1),tC(ceil(length(tC)/zoomDisp))])
end
sgtitle('Forecasted Trajectory using Physics Discrepancy Model')
legend('True','Augmented','Plato')
hold off
set(gcf,'position',[300,100,800,600])

%% Calculate RSME of dynamical error
% 
% clear RMSE_plato RMSE_aug perChange
% 
% for i = 1:nOrder
%     RMSE_plato(:,i) = sqrt(mean((xA(:,i) - xC(:,i)).^2));
%     RMSE_aug(:,i) = sqrt(mean((xA(:,i) - xB(:,i)).^2));
%     perChange(i) = (RMSE_plato(:,i)-RMSE_aug(:,i))/(RMSE_plato(:,i))*100;
% end
% 
% disp(['RMSE of Platonic dynamical error: ',num2str(RMSE_plato)])
% disp(['RMSE of Augemented dynamical error: ',num2str(RMSE_aug)])
% disp(['Percent decrease in RMSE: ',num2str(perChange),'%'])

return;