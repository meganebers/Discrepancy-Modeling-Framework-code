clear all, close all, clc
addpath('./DYNAMICS'); % system dynamics
addpath('./brunton_DMDbook/CODE/CH01_INTRO'); 
addpath('./util'); % other functions
addpath('./sparsedynamics - Brunton/utils'); 

%% Generage true dynamics (z) and measurement dynamics (y)

noise = 0.1;
dt = 0.1;
tlength = 60;
tspan = [0:dt:tlength];

system = 'Vanderpol';

eval(['discrepancyDynamics_',system])

T = tx;
D = delta;
X = x;
Y = y;
Z = z;

%% plot figure

figure,
plot(T,Y,'k','Linewidth',1),
hold on,
plot(T,Z,'r--','Linewidth',1),
hold on,
plot(T,X,'b--','Linewidth',1),
legend('Noisy Measurements','True System','Inadequate System'),
xlabel('Time'),
ylabel('Amplitude'),
title(''),
hold off


figure,
plot(Y(:,1),Y(:,2),'k','Linewidth',1),
hold on,
plot(Z(:,1),Z(:,2),'r--','Linewidth',1),
hold on,
plot(X(:,1),X(:,2),'b--','Linewidth',1),
legend('Noisy Measurements','True System','Inadequate System'),
xlabel('Time'),
ylabel('Amplitude'),
title(''),
hold off

