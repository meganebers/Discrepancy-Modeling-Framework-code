%% Generate data from ideal toy model

n = 3;
%addState = n;

%% lorenz system

sigma = 10;  % Lorenz's parameters (chaotic)
beta = 8/3;
rho = 28;
x0 = [-8; 8; 27];  % Initial condition
ep = 0.1; % discrepancy magnitude

N = length(tspan);
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));

[tx,x] = ode45(@(t,x) lorenz(t,x,sigma,beta,rho),tspan,x0,options); % ideal model data

figure,
for i = 1:n
    subplot(1,n,i)
    plot(tx,x(:,i),'k','LineWidth',1), hold on
end

%% Generate true dynamics (simulated) and measured states

[tz,z] = ode45(@(t,z) lorenz_discrep(t,z,sigma,beta,rho,g),tspan,x0,options); % true model

ty = tz;
y = z + noise*randn(size(z));
 
for i = 1:n
    subplot(1,n,i)
    plot(ty,y(:,i),'r--','LineWidth',1), hold on
end
xlabel('Time')
sgtitle('Measured Lorenz Dynamics')
hold off

figure, 
plot3(y(:,1),y(:,2),y(:,3),'r--','LineWidth',1), grid on,
hold on,
plot3(x(:,1),x(:,2),x(:,3),'k','LineWidth',1), grid on,
legend('Measured','Ideal','Location','Northeast')
sgtitle('Lorenz System')
xlabel('x'),ylabel('y')
hold off


%% Derivative of noisy measured dynamics; must smooth signal before numerical differentiation

%{
for i = 1:n
    y_butter(:,i) = butterfilt( y(:,i),cutoff,order, 'Low',ty,'on'  );
end

E = y_butter - x; % SS error
Eclean = z - x;

%dE(:,1) = gradient(E(:,1),dt);
%dE(:,2) = gradient(E(:,2),dt);

dy(:,1) = gradient(y_butter(:,1),dt);
dy(:,2) = gradient(y_butter(:,2),dt);

for i=1:length(x)
    dx(i,:) = vanderpol(0,y_butter(i,:),mu); % dx = f(x)
    dx_clean(i,:) = vanderpol(0,z(i,:),mu); % dx = f(x)

end

Ef = dy - dx; % discrepancy signal
Ef_clean = dz - dx_clean; % discrepancy signal

figure,
subplot(1,2,1),plot(ty,Ef,'Linewidth',[1]), grid on,
title('Identified Discrepancy (from measurements)')
subplot(1,2,2),plot(tz,Ef_clean,'Linewidth',[1]), grid on,
title('True Discrepancy (for comparison)')
%ylim([-10,10])
xlabel('time'), ylabel('amplitude')
sgtitle('Discrepancy Dynamics')

%}
%% compute Derivative

clear dy dx dxclean dz e ef efclean

if lowpass_filter == 1
    
    y_tmp = [];
    for i = 1:n
        y_tmp(:,i) = butterfilt( y(:,i),cutoff,4, 'Low',ty,'on'  );
    end
    
    y = y_tmp;
    
    e = y - x; % filtering
    
else
    
    e = y - x; % no filtering
    
end

for i = 1:n
    dy(:,i) = gradient(y(:,i),dt);
end

% Derivative of true dynamics dz = f(x) + g(x), (no noise) 
% (just for comparison!)
for i=1:length(z)
    dz(i,:) = lorenz_discrep(0,z(i,:),sigma,beta,rho,g); 
end

% Detrivative of ideal dynamics *I know what dx = f(x) is so I can plug my measurements directly into my model to find my derivatives)
% to be used in \delta\Phi = dy - f(x) = \Theta\Xi
for i=1:length(x)
    dx(i,:) = lorenz(0,y(i,:),sigma,beta,rho); % dx = f(x)
    dxclean(i,:) = lorenz(0,z(i,:),sigma,beta,rho); % dx = f(x)
end

% Subtract signals to isolate discrepancy 
ef = dy - dx;
efclean = dz - dxclean;
