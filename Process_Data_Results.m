%% Process data results from Generate_All_Data.m

% Written by Megan Ebers at the University of Washington 2021
% code to accompany paper submitted to SIAM-DS
% Discrepancy Modeling Framework: Learning missing physics, modeling 
% systematic residuals, and disambiguating between deterministic and random effects
% https://arxiv.org/abs/2203.05164

%% Make figure for showing Vanderpol Oscillator + discrepancy

clear all, close all, clc,
set(0,'DefaultFigureVisible','on')

system = 'Lorenz';
addpath('./util')
load([system,'_All_Results_2.mat'])

y = DMdata.IDphysics.NN.nonoise.y;
x = DMdata.IDphysics.NN.nonoise.x;

t = (1:length(y))*0.01;

switch system
    case  'Lorenz'
    figure, subplot(1,4,1), plot3(x(:,1), x(:,2), x(:,3), 'Linewidth',[2]), hold on,  plot3(y(:,1), y(:,2), y(:,3), 'Linewidth',[2]), grid off
    title('Phase Portrait'), view([35 10])
    case 'Vanderpol'
    figure, subplot(1,4,1), plot(x(:,1), x(:,2), 'Linewidth',[2]), hold on,plot(y(:,1), y(:,2), 'Linewidth',[2]), title('Phase Portrait'), grid off
end

subplot(1,4,[2 3]), plot(t,x(:,1), 'Linewidth',[4]), hold on, plot(t,y(:,1), 'Linewidth',[4]), grid on
xlim([t(1) t(end)]),title('Time Evolution'), %legend('$x$','$y_k$','Interpreter','Latex', 'FontSize',[9])
subplot(1,4,4), plot(t,y(:,1)-x(:,1), 'Linewidth',[2]), grid on,
title('Discrepancy')
sgtitle([system,' True vs. Approximate'])
set(gcf,'position',[100,300,1400,300],'color','w')
%saveas(gcf,['./Results/',system,'/',system,'_System_Dynamics.svg'])

%% Doubel check Forecasted Augemented dynamics are correct

methods = {'SINDy','DMD','GPR','NN'};
approach = {'IDphysics','errorModel'};
noiseLvl = {'mednoise'};

figure,
count = 1;
for j = 1:length(noiseLvl)
    for i = 1:length(approach)
        for k = 1:length(methods)
            
            eval(['dataset = DMdata.',approach{i},'.',methods{k},'.',noiseLvl{j},';']);
            
            xA = dataset.xA;
            xB = dataset.xB;
            xC = dataset.xC;
            
            subplot(1,length(methods)*length(approach)*length(noiseLvl),count), plot(xB), hold on, plot(xA)
            title([approach{i},' ',methods{k}])
            
            count = count + 1;
            
        end
    end
end

%% Update xB in DMD error Model for all noise levels

methods = {'DMD'};
approach = {'errorModel'};
noiseLvl = {'nonoise','lownoise','mednoise','highnoise'};

figure,
count = 1;
for j = 1:length(noiseLvl)
    for i = 1:length(approach)
        for k = 1:length(methods)
            
            eval(['dataset = DMdata.',approach{i},'.',methods{k},'.',noiseLvl{j},';']);
            
            xA = dataset.xA;
            xB = dataset.xB;
            xC = dataset.xC;
            
            %xB = xB(:,1:n) + xC;
            
            subplot(1,length(methods)*length(approach)*length(noiseLvl),count), plot(xB), hold on, plot(xA)
            title([approach{i},' ',methods{k}])
            
            count = count + 1;
            
        end
    end
end

%% Calculate RMSE for Training Data Reconstruction

methods = {'SINDy','DMD','GPR','NN'};
approach = {'IDphysics','errorModel'};
noiseLvl = {'nonoise','lownoise','mednoise','highnoise'};
vars = {'x','y','xR','time'};

figure,

count = 0;
dt = 0.01;
%RMSEwindow = 3;

for j = 1:length(noiseLvl)
    for i = 1:length(approach)
        for k = 1:length(methods)
            
            eval(['dataset = DMdata.',approach{i},'.',methods{k},'.',noiseLvl{j},';'])
            
            switch system
                
                case 'Vanderpol'
                    
                    [row,col]=size(dataset.xR);
                    [RMSE_plato, RMSE_aug, perChange] = Compare_RSME(dataset.y(1:row,:),dataset.xR(1:row,:),dataset.x(1:row,:),1)%,RMSEwindow/dt);
                    
                case 'Lorenz'
                    
                    [row,col]=size(dataset.xR);
                    [RMSE_plato, RMSE_aug, perChange] = Compare_RSME(dataset.y(1:row,:),dataset.xR(1:row,:),dataset.x(1:row,:),col)%,RMSEwindow/dt);
                    
            end
            
            residual.plato{k} = RMSE_plato;
            residual.aug{k} = RMSE_aug;
            residual.perChange{k} = perChange;
            
        end
        
        count = count + 1;
        
        subplot(length(noiseLvl),length(approach),count),
        X = categorical(methods);
        X = reordercats(X,methods);
        Y = [median(residual.perChange{1}); median(residual.perChange{2}); median(residual.perChange{3}); median(residual.perChange{4})];
        bar(X,Y)
        ylim([-100, 100])
        title(noiseLvl{j})
        ylabel('RMSE Decrease (%)'), hold on,
        set(gcf,'position',[1000,0,500,800],'color','w')
        
    end
end

sgtitle([system,' Reconstruction'])

%% Calculate RMSE for Test Data Forecasting

methods = {'SINDy','DMD','GPR','NN'};
approach = {'IDphysics','errorModel'};
noiseLvl = {'nonoise','lownoise','mednoise','highnoise'};
vars = {'xA','xB','xC','time'};

figure,

count = 0;
dt = 0.01;
%RMSEwindow = 3;

for j = 1:length(noiseLvl)
    for i = 1:length(approach)
        for k = 1:length(methods)
            
            eval(['dataset = DMdata.',approach{i},'.',methods{k},'.',noiseLvl{j},';'])
            
            switch system
                
                case 'Vanderpol'
                    
                    [row,col]=size(dataset.xA);
                    [RMSE_plato, RMSE_aug, perChange] = Compare_RSME(dataset.xA,dataset.xB,dataset.xC,1)%,RMSEwindow/dt);
                    
                case 'Lorenz'
                    
                    [row,col]=size(dataset.xA);
                    [RMSE_plato, RMSE_aug, perChange] = Compare_RSME(dataset.xA,dataset.xB,dataset.xC,col)%,RMSEwindow/dt);
                    
            end
            
            residual.plato{k} = RMSE_plato;
            residual.aug{k} = RMSE_aug;
            residual.perChange{k} = perChange;
            
        end
        
        count = count + 1;
        
        subplot(length(noiseLvl),length(approach),count),
        X = categorical(methods);
        X = reordercats(X,methods);
        Y = [median(residual.perChange{1}); median(residual.perChange{2}); median(residual.perChange{3}); median(residual.perChange{4})];
        bar(X,Y)
        ylim([-100, 100])
        title(noiseLvl{j})
        ylabel('RMSE Decrease (%)'), hold on,
        set(gcf,'position',[1000,0,500,800],'color','w')
        
    end
    
end

sgtitle([system,' Forecasting'])

return;

%% Time Cost for all

clear all, close all, clc,
set(0,'DefaultFigureVisible','on')

system = 'Lorenz';
addpath('./util')
load([system,'_All_Results_3.mat'])

methods = {'SINDy','DMD','GPR','NN'};
approach = {'IDphysics','errorModel'};
noiseLvl = {'nonoise','lownoise','mednoise','highnoise'};
vars = {'xA','xB','xC','time'};

figure,

count = 0;
dt = 0.01;

for j = 1:length(noiseLvl)
    for i = 1:length(approach)
        for k = 1:length(methods)
            
            eval(['dataset = DMdata.',approach{i},'.',methods{k},'.',noiseLvl{j},';'])
            
            residual.time{k} = dataset.time;
            
        end
        
        count = count + 1;
        
        subplot(length(noiseLvl),length(approach),count),
        X = categorical(methods);
        X = reordercats(X,methods);
        Y = [residual.time{1}; residual.time{2}; residual.time{3}; residual.time{4}];
        bar(X,Y)
        %ylim([0,  500])
        title(noiseLvl{j})
        %ylabel('Computational Cost (seconds)'), hold on,
        
    end
end

sgtitle([system,' Computational Cost'])
set(gcf,'position',[1000,0,500,800],'color','w'), hold on,

%saveas(gcf,['./Results/',system,'/',system,'_Comp_Cost.svg'])

%% Plot variation of RMSE percent change over testSet time

clear all, close all, clc,
set(0,'DefaultFigureVisible','on')

system = 'Lorenz';
addpath('./util')
load([system,'_All_Results_2.mat'])

xA = DMdata.IDphysics.NN.nonoise.xA;
xB = DMdata.IDphysics.NN.nonoise.xB;
xC = DMdata.IDphysics.NN.nonoise.xC;

dt = 0.01;
n = 3;

for j = 1:1:20
    [RMSE_plato, RMSE_aug, perChange(j,:)] = Compare_RSME(xA,xB,xC,n,j/dt);
end
figure, plot(1:1:20,perChange,'Linewidth',[2]), hold on,
yline([0]),grid on,
legend('x','y','z')
xlabel('Window of Forecasting (seconds)')
ylabel('% Change in RMSE (positive = better performance)')
title('Effectiveness of Discrepancy Model Augmentation Over Longer Forecasting Windows')
xlim([1,20])
set(gcf,'position',[0,0,1000,400],'color','w')

tA = (1:length(xA))*dt;

figure,
for i = 1:n
    subplot(1,n,i)
    plot(tA,xA(:,i),'k','Linewidth',[2]), hold on % Truth
    %plot(tA,xC(:,i),'b--','Linewidth',[2]), hold on % Plato
    plot(tA,xB(:,i),'r--','Linewidth',[2]), hold on % Augmented
    xlim([tA(1),tA(roundn(length(tA), 2)/2)])
    grid on,
end
sgtitle('Forecasted Trajectory using Discrepancy Model')
%legend('True','Augmented')

set(gcf,'position',[100,300,1200,400],'color','w')

%% Visualize RMSE

figure,
subplot(1,2,1),
X = categorical({'SINDy','DMD','GP','NN'});
X = reordercats(X,{'SINDy','DMD','GP','NN'});
Y1 = [median(perChange_TOTAL.SINDy_IDphysics); median(perChange_TOTAL.DMD_IDphysics); median(perChange_TOTAL.GP_IDphysics); median(perChange_TOTAL.NN_IDphysics)];
bar(X,Y1)
ylim([-100, max(Y1)])
title('Percent RMSE Decrease for Learning Physics')
ylabel('RMSE')

subplot(1,2,2),
X = categorical({'SINDy','DMD','GP','NN'});
X = reordercats(X,{'SINDy','DMD','GP','NN'});
Y2 = [perChange_TOTAL.SINDy_errorModel,perChange_TOTAL.DMD_errorModel,perChange_TOTAL.GP_errorModel,perChange_TOTAL.NN_errorModel];
bar(X,Y2)
ylim([-100, max(Y2)])
title('Percent RMSE Decrease for Modeling Error')
ylabel('RMSE')

sgtitle(['Comparing Discrepancy Model Performance for the ',system,' system with ',num2str(noise*100),'% Noise for ',num2str(RMSEwindow),' seconds'])
set(gcf,'position',[100,300,1000,600])

%% Compare each method for each approach - NO NOISE ONLY

clear all, close all, clc,
set(0,'DefaultFigureVisible','on')

system = 'Vanderpol';
n = 1;
addpath('./util')
load([system,'_All_Results_2.mat'])

methods = {'SINDy','DMD','GPR','NN'};
approach = {'IDphysics','errorModel'};
noiseLvl = {'highnoise'};

for j = 1:length(noiseLvl)
    for i = 1:length(approach)
        for k = 1:length(methods)
            
            eval(['dataset = DMdata.',approach{i},'.',methods{k},'.',noiseLvl{j},';'])
            
            x = dataset.x(:,n);
            y = dataset.y(:,n);
            xR = dataset.xR(:,n);
            xA = dataset.xA(:,n);
            xB = dataset.xB(:,n);
            xC = dataset.xC(:,n);
            
            
            dt = 0.01;
            tR = (1:length(x))*dt;
            t_tmp = (1:length(xA))*dt;
            tF = ones(1,length(xA))*length(x)*dt + t_tmp;
            
            % Plot
            
            figure,
            % Reconstruction
            plot(tR,y,'r','Linewidth',[2]), hold on,
            plot(tR,x,'b','Linewidth',[2]), hold on,
            plot(tR,xR,'k--','Linewidth',[2]), hold on,
            % Forecasting
            plot(tF,xA,'r','Linewidth',[2]), hold on,
            plot(tF,xC,'b','Linewidth',[2]), hold on,
            plot(tF,xB,'k--','Linewidth',[2]), hold on,
            xline(tF(1),'k','Linewidth',[1])
            
            xlim([0 tF(end)])
            grid on
            legend('y','x','$\tilde{x}$','Interpreter','Latex')
            title(['Waveform of ',system,' using ',approach{i},' - ',methods{k},' (',noiseLvl{j},')'])
            xlabel('Time (s)')
            ylabel('Amplitude')
            set(gcf,'position',[100,300,1200,400])
            
            saveas(gcf,['./Results/Waveforms/',system,'/',noiseLvl{j},'/',system,'_',approach{i},'_',methods{k},'_',noiseLvl{j},'.svg'])
            
        end
    end
end

%% Compare each method for each approach - COMPARING ERROR

clear all, close all, clc,
set(0,'DefaultFigureVisible','on')

system = 'Lorenz';
n = 1;
addpath('./util')
load([system,'_All_Results_2.mat'])

methods = {'SINDy','DMD','GPR','NN'};
approach = {'IDphysics','errorModel'};
noiseLvl = {'highnoise'};

for j = 1:length(noiseLvl)
    for i = 1:length(approach)
        for k = 1:length(methods)
            
            eval(['dataset = DMdata.',approach{i},'.',methods{k},'.',noiseLvl{j},';'])
            
            x = dataset.x(:,n);
            y = dataset.y(:,n);
            xR = dataset.xR(:,n);
            xA = dataset.xA(:,n);
            xB = dataset.xB(:,n);
            xC = dataset.xC(:,n);
            
            
            dt = 0.01;
            tR = (1:length(x))*dt;
            t_tmp = (1:length(xA))*dt;
            tF = ones(1,length(xA))*length(x)*dt + t_tmp;
            
            % Plot
            
            figure,
            % Reconstruction
            plot(tR,y-x,'b','Linewidth',[2]), hold on,
            plot(tR,y-xR,'k--','Linewidth',[2]), hold on,
            % Forecasting
            plot(tF,xA-xC,'b','Linewidth',[2]), hold on,
            plot(tF,xA-xB,'k--','Linewidth',[2]), hold on,
            xline(tF(1),'k','Linewidth',[1])
            
            xlim([0 tF(end)])
            grid on
            legend('y-x','$y-\tilde{x}$','Interpreter','Latex')
            title(['Waveform of ',system,' using ',approach{i},' - ',methods{k},' (',noiseLvl{j},')'])
            xlabel('Time (s)')
            ylabel('Amplitude')
            set(gcf,'position',[100,300,1200,400])
            
            saveas(gcf,['./Results/Waveforms/',system,'/error/',noiseLvl{j},'/',system,'_',approach{i},'_',methods{k},'_',noiseLvl{j},'.svg'])
            
        end
    end
end

%% Make waveform of y and x with phase portrait too
clear all, close all, clc,
set(0,'DefaultFigureVisible','on')

system = 'Vanderpol';
n = 2;
addpath('./util')
load([system,'_All_Results_2.mat'])

eval(['dataset = DMdata.IDphysics.SINDy.nonoise;'])
x = dataset.x(:,n);
y = dataset.y(:,n);
xR = dataset.xR(:,n);
xA = dataset.xA(:,n);
xB = dataset.xB(:,n);
xC = dataset.xC(:,n);

