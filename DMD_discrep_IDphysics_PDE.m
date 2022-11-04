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

clear all, close all, clc
addpath('./DYNAMICS'); % system dynamics
addpath('./util'); % other functions
addpath('./optdmd-master/src'); % optDMD
addpath('./ODE_Solvers');
addpath('./sparsedynamics - Brunton/utils');

%% Generage true dynamics (z) and measurement dynamics (y)
%
noise = 0.0000;
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
%}

system = 'Burgers';

eval(['discrepancyDynamics_',system])

T = tx.';
E = e.'; % state space error
Ef = ef.'; % dynamical model error
X = real(x).';
Y = real(y).';
Z = real(z).';

trainSet.T = T(1:round(PD*length(T))); 
testSet.T = T(round(PD*length(T)):end);

trainSet.E = E(1:round(PD*length(E)),:);
testSet.E = E(round(PD*length(E)):end,:);

trainSet.Ef = Ef(1:round(PD*length(Ef)),:);
testSet.Ef = Ef(round(PD*length(Ef)):end,:);

trainSet.X = X(1:round(PD*length(X)),:);
testSet.X = X(round(PD*length(X)):end,:);

trainSet.Y = Y(1:round(PD*length(Y)),:);
testSet.Y = Y(round(PD*length(Y)):end,:);

trainSet.Z = Z(1:round(PD*length(Z)),:);
testSet.Z = Z(round(PD*length(Z)):end,:);

%% Perform optDMD

r = 15;

% linear constraints; meant to constrain eigs to neg real
lbc = [-Inf*ones(r); -Inf*ones(r)];
ubc = [zeros(r); Inf*ones(r)];

copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);

ts = trainSet.T;
imode = 2;
[DMD.Phi,DMD.lambda,DMD.b] = optdmd(trainSet.Ef.',ts,r,imode,[],[],[],copts);
DMD.omega = log(DMD.lambda)/dt; % discrete to continuous time eigs
disp('DMD computed')

%% Reconstruction of the discrepancy

%D_recon = w*diag(b)*exp(e*t')
D1 = DMD.Phi*diag(DMD.b)*exp(DMD.lambda*ts);
%relerr_r = norm(D1-Y_H,'fro')/norm(Y_H,'fro')
%relerr_r_clean = norm(D1-xclean,'fro')/norm(xclean,'fro');

figure,
subplot(1,2,1), waterfall(xgrid,ts,real(trainSet.Ef)), 
view(0,90), set(gca,'Fontsize',[12]), hold on,
subplot(1,2,2), waterfall(xgrid,ts,real(D1.')),
view(0,90), set(gca,'Fontsize',[12]), hold on,
grid on,
colormap(flipud(jet));
sgtitle('DMD Reconstruction of First Time Delay')
set(gcf,'position',[300,100,1000,400])
hold off

%% Reconstruction of Training Data

tspan = trainSet.T;

% setup
eps=0.1;
L=16; n=256;
x2=linspace(-L/2,L/2,n+1); xgrid=x2(1:n); k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
k2=fftshift(k);

% initial data
%u=trainSet.Y(end,:).';
u=exp(-(xgrid+2).^2).';
ut=fft(u);

xR = ut;
for i = 1:length(tspan)-1
    tt = [tspan(i) tspan(i+1)];
    xR_tmp = ode5('burgers_approx_DMD',tt,ut,k,eps,DMD,dt,tt(1));
    ut = xR_tmp(2,:).';
    xR = [xR ut];
end
xR_tmp = xR.';
tR =tspan;

xR = [];
for j=1:length(tspan)
    xR(:,j)=ifft( xR_tmp(j,1:n).' );
end

%% Plotting Reconstruction Results

tt = tR;
data_state = {real(trainSet.Y) real(trainSet.X) real(xR.')};
data_error = {real(trainSet.Y-trainSet.Y) real(trainSet.Y-trainSet.X) real(trainSet.Y-xR.')};

% Reconstruction states

figure,
ax(1) = subplot(1,4,1);
waterfall(xgrid,tt,data_state{1});
ylabel('Time')
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2])
title('True'), hold on,
colormap(ax(1), b2r( min(min([data_state{:}])),max(max([data_state{:}])) ) );
colorbar;

ax(2) = subplot(1,4,2);
waterfall(xgrid,tt,data_state{2});
ylabel('Time')
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2])
title('Plato'), hold on,
colormap(ax(2), b2r( min(min([data_state{:}])),max(max([data_state{:}])) ) );
colorbar;

ax(3) = subplot(1,4,3);
waterfall(xgrid,tt,data_state{3});
ylabel('Time')
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2])
title('Augmented'), hold on,
colormap(ax(3), b2r( min(min([data_state{:}])),max(max([data_state{:}])) ) );
colorbar;

ax(4) = subplot(1,4,4);
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2])
colormap(ax(4), b2r( min(min([data_state{:}])),max(max([data_state{:}])) ) ); 
colorbar;

sgtitle('Reconstruction States using Discrepancy Model')
set(gcf,'position',[100,100,1500,300])

%saveas(gcf,['./Results/Burgers/Burgers_Reconstruction_States.fig'])


%% Reconstruction error 

figure,
ax(1) = subplot(1,4,1);
waterfall(xgrid,tt,data_error{1}); 
ylabel('Time')
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2])
title('True - Resiudal Error'), hold on,
colormap(ax(1), b2r( -eps, eps ))
colorbar;

ax(2) = subplot(1,4,2);
waterfall(xgrid,tt,data_error{2}); 
ylabel('Time')
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2])
title('Plato - Resiudal Error'), hold on,
colormap(ax(2), b2r( min(min([data_error{2}])),max(max([data_error{2}])) )),
colorbar;

ax(3) = subplot(1,4,3);
waterfall(xgrid,tt,data_error{3});
ylabel('Time')
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2])
title('Augmented - Resiudal Error'), hold on,
colormap(ax(3), b2r( min(min([data_error{3}])),max(max([data_error{3}])) ))
colorbar;

ax(4) = subplot(1,4,4);
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2])
colormap(ax(4), b2r( min(min([data_error{:}])),max(max([data_error{:}])) )); 
colorbar;

sgtitle('Reconstruction Error using Discrepancy Model')
set(gcf,'position',[100,100,1500,300])

%saveas(gcf,['./Results/Burgers/Burgers_Reconstruction_Error.fig'])

%% combine DMD model and toy system for forecasting

tspan = testSet.T;
% ***** MUST START AT 0 due to discrepancy error == 0 at initialization 
%       regardless of IC, but DMD is dependent on time series and needs
%       to be initialized at t = 0;

% setup
eps=0.1;
L=16; n=256;
x2=linspace(-L/2,L/2,n+1); xgrid=x2(1:n); k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
k2=fftshift(k);

% initial data
%u=trainSet.Y(end,:).';
u=exp(-(xgrid+2).^2).';
ut=fft(u);

[tA,xA_tmp]=ode45('burgers_rhs_discrep',tspan,ut,[],k,eps); % true model in test domain

if lowpass_filter == 1
    for i = 1:n
        xA_tmp(:,i) = butterfilt( xA_tmp(:,i),cutoff,4, 'Low',ty,'off'  ); % filtering
    end
    xA_tmp = xA_tmp;
else
    xA_tmp = xA_tmp; % no filtering
end

xB = ut;
for i = 1:length(tspan)-1
    tt = [tspan(i) tspan(i+1)];
    xB_tmp = ode5('burgers_approx_DMD',tt,ut,k,eps,DMD,dt,tt(1));
    ut = xB_tmp(2,:).';
    xB = [xB ut];
end
xB_tmp = xB.';
tB =tspan;

[tC,xC_tmp]=ode45('burgers_rhs',tspan,ut,[],k,eps); % ideal model (inadequate)

xA = [];
xB = [];
xC = [];
for j=1:length(tspan)
    xA(:,j)=ifft( xA_tmp(j,1:n).' );
    xB(:,j)=ifft( xB_tmp(j,1:n).' );
    xC(:,j)=ifft( xC_tmp(j,1:n).' );
end

%% Plotting Forecasting Results
 
tt = tA;
data_state = {real(xA.') real(xC.') real(xB.')};
data_error = {real(xA.'-xA.') real(xA.'-xC.') real(xA.'-xB.')};

% Forecasting states

figure,
ax(1) = subplot(1,4,1);
waterfall(xgrid,tt,data_state{1});
ylabel('Time')
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2])
title('True'), hold on,
colormap(ax(1), b2r( min(min([data_state{:}])),max(max([data_state{:}])) ) );
colorbar;

ax(2) = subplot(1,4,2);
waterfall(xgrid,tt,data_state{2});
ylabel('Time')
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2])
title('Plato'), hold on,
colormap(ax(2), b2r( min(min([data_state{:}])),max(max([data_state{:}])) ) );
colorbar;

ax(3) = subplot(1,4,3);
waterfall(xgrid,tt,data_state{3});
ylabel('Time')
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2])
title('Augmented'), hold on,
colormap(ax(3), b2r( min(min([data_state{:}])),max(max([data_state{:}])) ) );
colorbar;

ax(4) = subplot(1,4,4);
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2])
colormap(ax(4), b2r( min(min([data_state{:}])),max(max([data_state{:}])) ) ); 
colorbar;

sgtitle('Forecasted States using Discrepancy Model')
set(gcf,'position',[100,100,1500,300])

%saveas(gcf,['./Results/Burgers/Burgers_Forecasted_States.fig'])


% Forecasting error 

figure,
ax(1) = subplot(1,4,1);
waterfall(xgrid,tt,data_error{1}); 
ylabel('Time')
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2])
title('True - Resiudal Error'), hold on,
colormap(ax(1), b2r( -eps, eps ))
colorbar;

ax(2) = subplot(1,4,2);
waterfall(xgrid,tt,data_error{2}); 
ylabel('Time')
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2])
title('Plato - Resiudal Error'), hold on,
colormap(ax(2), b2r( min(min([data_error{:}])),max(max([data_error{:}])) )),
colorbar;

ax(3) = subplot(1,4,3);
waterfall(xgrid,tt,data_error{3});
ylabel('Time')
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2])
title('Augmented - Resiudal Error'), hold on,
colormap(ax(3), b2r( min(min([data_error{:}])),max(max([data_error{:}])) ))
colorbar;

ax(4) = subplot(1,4,4);
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2])
colormap(ax(4), b2r( min(min([data_error{:}])),max(max([data_error{:}])) )); 
colorbar;

sgtitle('Forecasted Error using Discrepancy Model')
set(gcf,'position',[100,100,1500,300])

%saveas(gcf,['./Results/Burgers/Burgers_Forecasting_Error.fig'])

%% Calculate RSME of dynamical error

clear RMSE_plato RMSE_aug perChange

for i = 1:length(testSet.T)
    RMSE_plato(:,i) = sqrt( mean( real( xA(:,i) ) - real( xC(:,i) ) ).^2);
    RMSE_aug(:,i) = sqrt( mean( real( xA(:,i) ) - real( xB(:,i) ) ).^2);
    perChange(i) = ( RMSE_plato(:,i) - RMSE_aug(:,i) ) / ( RMSE_plato(:,i) ) * 100;
end

disp(['RMSE of Platonic dynamical error: ',num2str(RMSE_plato)])
disp(['RMSE of Augemented dynamical error: ',num2str(RMSE_aug)])
disp(['Percent decrease in RMSE: ',num2str(perChange),'%'])

figure, plot(perChange,'Linewidth',[2]), 
xlabel('Time'), ylabel('Percent Decreased Error'), sgtitle('Change in Percent Decreased Error over Time')
