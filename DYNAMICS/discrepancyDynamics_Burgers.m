%% Add path
addpath('./DYNAMICS/PDE'); % system dynamics

%% Burgers' equation
%  u_t + u*u_x - eps*u_xx =0 

n = 1;

% setup
eps=0.1; 
L=16; n=256; 
x2=linspace(-L/2,L/2,n+1); xgrid=x2(1:n); k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
k2=fftshift(k);

% initial data
u=exp(-(xgrid+2).^2).';
ut=fft(u); 

[tx,utsol]=ode45('burgers_rhs',tspan,ut,[],k,eps);
  
x = [];

for j=1:length(tspan)
    x(:,j)=ifft( utsol(j,1:n).' );
end

figure(1), 
subplot(1,3,1),
C = gradient(real(x.'));
waterfall(xgrid,tx,real(x.'),C);
ylabel('Time')
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2],'Ylim',[0 tspan(end)],'Ytick',[0 10])
title('Ideal')

%% Generate true dynamics (simulated) and measured states

[tz,ztsol]=ode45('burgers_rhs_discrep',tspan,ut,[],k,eps);

ty = tz;

z = [];
for j=1:length(tspan)
    z(:,j)=ifft( ztsol(j,1:n).' );
end

y = z + noise*randn(size(z));

figure(1), 
subplot(1,3,2),
C = gradient(real(y.'));
waterfall(xgrid,ty,real(y.'),C); 
view(0,90), set(gca,'Fontsize',[12])
xlabel('Space'), 
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2],'Ylim',[0 tspan(end)],'Ytick',[0 10])
title('Measured')

%% compute Derivative

clear dy dx dxclean dz E Ef Efclean

if lowpass_filter == 1
    
    for i = 1:length(tspan)
        y_tmp(:,i) = butterfilt( real(y(:,i)),cutoff,4, 'Low',ty,'off'  );
    end
    y = y_tmp;
    
     e = real(y) - real(x); % filtering
    
else
    
    e = (real(y) - real(x)); % no filtering
    
end

figure(1), 
subplot(1,3,3),
C = gradient(real(e.'));
waterfall(xgrid,ty,real(e.'),C); 
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2],'Ylim',[0 tspan(end)],'Ytick',[0 10])
title('Discrepancy')
set(gcf,'position',[300,100,1000,400])
hold off

[px,py] = gradient(real(y),dt);
dy = px;

% Derivative of true dynamics dz = f(x) + g(x), (no noise) 
% (just for comparison!)
zt=fft(real(z));
for j=1:length(z)
    dzt(:,j) = burgers_rhs_discrep(0,zt(:,j),[],k,eps); 
end

for j=1:length(tspan)
    dz(:,j)=ifft( dzt(:,j) );
end

% Detrivative of ideal dynamics *I know what dx = f(x) is so I can plug my measurements directly into my model to find my derivatives)
% to be used in \delta\Phi = dy - f(x) = \Theta\Xi
zt=fft(real(z));
yt=fft(real(y));
for j=1:length(z)
    dxt(:,j) = burgers_rhs(0,yt(:,j),[],k,eps);
    dxtclean(:,j) = burgers_rhs(0,zt(:,j),[],k,eps);
end

for j=1:length(tspan)
    dx(:,j)=ifft( dxt(:,j) );
    dxclean(:,j)=ifft( dxtclean(:,j) );
end

% Subtract signals to isolate discrepancy 
ef = (real(dy) - real(dx));
efclean = (real(dz) - real(dxclean));

figure(2), 
subplot(1,2,1),
%C = gradient(real(E.'));
waterfall(xgrid,ty,real(e.')); 
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2],'Ylim',[0 tspan(end)],'Ytick',[0 10])
title('State Space Error'),
hold on,

subplot(1,2,2),
%C = gradient(real(Ef.'));
waterfall(xgrid,ty,real(ef.')); 
view(0,90), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2],'Ylim',[0 tspan(end)],'Ytick',[0 10])
title('Dynamical Error'),
set(gcf,'position',[300,100,1000,400])
hold off