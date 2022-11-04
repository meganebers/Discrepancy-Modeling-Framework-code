%% Generate all data

% Written by Megan Ebers at the University of Washington 2021
% code to accompany paper submitted to SIAM-DS
% Discrepancy Modeling Framework: Learning missing physics, modeling 
% systematic residuals, and disambiguating between deterministic and random effects
% https://arxiv.org/abs/2203.05164

%% Setup parameters

clear all, close all, clc,
set(0,'DefaultFigureVisible','off')

param.dt = 0.01;
param.tlength = 50;
param.tspan = [0:param.dt:param.tlength];
param.cutoff = 4;
param.order = 4;
param.PD = 0.60 ;
param.g = '0.01*y(1).*y(1).*y(1)'; % epsilon discrepancy

%{
% Vanderpol parameters
param.system = 'Vanderpol'; 
param.r = 15; % = rank
param.p = 100; % p = time delays, 
%}

% Lorenz parameters
param.system = 'Lorenz';  
param.r = 15; % = rank
param.p = 100; % p = time delays, 

%% SINDy Error Model

dt = param.dt;
tlength = param.tlength;
tspan = param.tspan;
cutoff = param.cutoff;
order = param.order;
PD = param.PD;
g = param.g;
system = param.system;
r = param.r;
p = param.p; 
    
noise = 0;
polyorder = 3; % polynomial order
usesine = 0; % 1 = yes, 0 = no
lambda = 0.046; % lambda is our sparsification knob
lowpass_filter = 0;

tStart = tic;
eval('SINDy_discrep_errorModel')
tEnd = toc(tStart);
    
DMdata.errorModel.SINDy.nonoise.xA = xA;
DMdata.errorModel.SINDy.nonoise.xB = xB;
DMdata.errorModel.SINDy.nonoise.xC = xC;
DMdata.errorModel.SINDy.nonoise.x = trainSet.X;
DMdata.errorModel.SINDy.nonoise.y = trainSet.Y;
DMdata.errorModel.SINDy.nonoise.xR = xR;
DMdata.errorModel.SINDy.nonoise.time = tEnd;

for noise = [0.001 0.01 0.1]
    
    clearvars -except param DMdata noise
    dt = param.dt;
    tlength = param.tlength;
    tspan = param.tspan;
    cutoff = param.cutoff;
    order = param.order;
    PD = param.PD;
    g = param.g;
    system = param.system;
    r = param.r;
    p = param.p;
    
    polyorder = 3; % polynomial order
    usesine = 0; % 1 = yes, 0 = no
    lambda = 0.046; % lambda is our sparsification knob
    lowpass_filter = 1;
    
    tStart = tic;
    eval('SINDy_discrep_errorModel')
    tEnd = toc(tStart);
    
    if noise == 0.001
        DMdata.errorModel.SINDy.lownoise.xA = xA;
        DMdata.errorModel.SINDy.lownoise.xB = xB;
        DMdata.errorModel.SINDy.lownoise.xC = xC;
        DMdata.errorModel.SINDy.lownoise.x = trainSet.X;
        DMdata.errorModel.SINDy.lownoise.y = trainSet.Y;
        DMdata.errorModel.SINDy.lownoise.xR = xR;
        DMdata.errorModel.SINDy.lownoise.time = tEnd;
    elseif noise == 0.01
        DMdata.errorModel.SINDy.mednoise.xA = xA;
        DMdata.errorModel.SINDy.mednoise.xB = xB;
        DMdata.errorModel.SINDy.mednoise.xC = xC;
        DMdata.errorModel.SINDy.mednoise.x = trainSet.X;
        DMdata.errorModel.SINDy.mednoise.y = trainSet.Y;
        DMdata.errorModel.SINDy.mednoise.xR = xR;
        DMdata.errorModel.SINDy.mednoise.time = tEnd;
    elseif noise == 0.1
        DMdata.errorModel.SINDy.highnoise.xA = xA;
        DMdata.errorModel.SINDy.highnoise.xB = xB;
        DMdata.errorModel.SINDy.highnoise.xC = xC;
        DMdata.errorModel.SINDy.highnoise.x = trainSet.X;
        DMdata.errorModel.SINDy.highnoise.y = trainSet.Y;
        DMdata.errorModel.SINDy.highnoise.xR = xR;
        DMdata.errorModel.SINDy.highnoise.time = tEnd;
    else
        disp('Please specifiy a noise level')
        return;
    end
    
end

disp('SINDy error model is done')

%% SINDy Learn Physics

clearvars -except param DMdata noise
dt = param.dt;
tlength = param.tlength;
tspan = param.tspan;
cutoff = param.cutoff;
order = param.order;
PD = param.PD;
g = param.g;
system = param.system;
r = param.r;
p = param.p;

noise = 0;
polyorder = 3; % polynomial order
usesine = 0; % 1 = yes, 0 = no
lambda = 0.0095; % lambda is our sparsification knob
lowpass_filter = 0;

tStart = tic;
eval('SINDy_discrep_IDphysics')
tEnd = toc(tStart);
    
DMdata.IDphysics.SINDy.nonoise.xA = xA;
DMdata.IDphysics.SINDy.nonoise.xB = xB;
DMdata.IDphysics.SINDy.nonoise.xC = xC;
DMdata.IDphysics.SINDy.nonoise.x = trainSet.X;
DMdata.IDphysics.SINDy.nonoise.y = trainSet.Y;
DMdata.IDphysics.SINDy.nonoise.xR = xR;
DMdata.IDphysics.SINDy.nonoise.time = tEnd;

for noise = [0.001 0.01 0.1]
    
    clearvars -except param DMdata noise
    dt = param.dt;
    tlength = param.tlength;
    tspan = param.tspan;
    cutoff = param.cutoff;
    order = param.order;
    PD = param.PD;
    g = param.g;
    system = param.system;
    r = param.r;
    p = param.p;
    
    polyorder = 3; % polynomial order
    usesine = 0; % 1 = yes, 0 = no
    lambda = 0.0095; % lambda is our sparsification knob
    lowpass_filter = 1;
    
    tStart = tic;
    eval('SINDy_discrep_IDphysics')
    tEnd = toc(tStart);
    
    if noise == 0.001
        DMdata.IDphysics.SINDy.lownoise.xA = xA;
        DMdata.IDphysics.SINDy.lownoise.xB = xB;
        DMdata.IDphysics.SINDy.lownoise.xC = xC;
        DMdata.IDphysics.SINDy.lownoise.x = trainSet.X;
        DMdata.IDphysics.SINDy.lownoise.y = trainSet.Y;
        DMdata.IDphysics.SINDy.lownoise.xR = xR;
        DMdata.IDphysics.SINDy.lownoise.time = tEnd;
    elseif noise == 0.01
        DMdata.IDphysics.SINDy.mednoise.xA = xA;
        DMdata.IDphysics.SINDy.mednoise.xB = xB;
        DMdata.IDphysics.SINDy.mednoise.xC = xC;
        DMdata.IDphysics.SINDy.mednoise.x = trainSet.X;
        DMdata.IDphysics.SINDy.mednoise.y = trainSet.Y;
        DMdata.IDphysics.SINDy.mednoise.xR = xR;
        DMdata.IDphysics.SINDy.mednoise.time = tEnd;
    elseif noise == 0.1
        DMdata.IDphysics.SINDy.highnoise.xA = xA;
        DMdata.IDphysics.SINDy.highnoise.xB = xB;
        DMdata.IDphysics.SINDy.highnoise.xC = xC;
        DMdata.IDphysics.SINDy.highnoise.x = trainSet.X;
        DMdata.IDphysics.SINDy.highnoise.y = trainSet.Y;
        DMdata.IDphysics.SINDy.highnoise.xR = xR;
        DMdata.IDphysics.SINDy.highnoise.time = tEnd;
    else
        disp('Please specifiy a noise level')
        return;
    end
    
end

disp('SINDy learn physics is done')

%% DMD Model Error
 
clearvars -except param DMdata noise
dt = param.dt;
tlength = param.tlength;
tspan = param.tspan;
cutoff = param.cutoff;
order = param.order;
PD = param.PD;
g = param.g;
system = param.system;
r = param.r;
p = param.p;

noise = 0;
lowpass_filter = 0;

tStart = tic;
eval('DMD_discrep_errorModel')
tEnd = toc(tStart);

DMdata.errorModel.DMD.nonoise.xA = xA;
DMdata.errorModel.DMD.nonoise.xB = xB;
DMdata.errorModel.DMD.nonoise.xC = xC;
DMdata.errorModel.DMD.nonoise.x = trainSet.X;
DMdata.errorModel.DMD.nonoise.y = trainSet.Y;
DMdata.errorModel.DMD.nonoise.xR = xR;
DMdata.errorModel.DMD.nonoise.time = tEnd;

for noise = [0.001 0.01 0.1]
    
    clearvars -except param DMdata noise
    dt = param.dt;
    tlength = param.tlength;
    tspan = param.tspan;
    cutoff = param.cutoff;
    order = param.order;
    PD = param.PD;
    g = param.g;
    system = param.system;
    r = param.r;
    p = param.p;
        
    lowpass_filter = 1;
    
    tStart = tic;
    eval('DMD_discrep_errorModel')
    tEnd = toc(tStart);
    
    if noise == 0.001
        DMdata.errorModel.DMD.lownoise.xA = xA;
        DMdata.errorModel.DMD.lownoise.xB = xB;
        DMdata.errorModel.DMD.lownoise.xC = xC;
        DMdata.errorModel.DMD.lownoise.x = trainSet.X;
        DMdata.errorModel.DMD.lownoise.y = trainSet.Y;
        DMdata.errorModel.DMD.lownoise.xR = xR;
        DMdata.errorModel.DMD.lownoise.time = tEnd;
    elseif noise == 0.01
        DMdata.errorModel.DMD.mednoise.xA = xA;
        DMdata.errorModel.DMD.mednoise.xB = xB;
        DMdata.errorModel.DMD.mednoise.xC = xC;
        DMdata.errorModel.DMD.mednoise.x = trainSet.X;
        DMdata.errorModel.DMD.mednoise.y = trainSet.Y;
        DMdata.errorModel.DMD.mednoise.xR = xR;
        DMdata.errorModel.DMD.mednoise.time = tEnd;
    elseif noise == 0.1
        DMdata.errorModel.DMD.highnoise.xA = xA;
        DMdata.errorModel.DMD.highnoise.xB = xB;
        DMdata.errorModel.DMD.highnoise.xC = xC;
        DMdata.errorModel.DMD.highnoise.x = trainSet.X;
        DMdata.errorModel.DMD.highnoise.y = trainSet.Y;
        DMdata.errorModel.DMD.highnoise.xR = xR;
        DMdata.errorModel.DMD.highnoise.time = tEnd;
    else
        disp('Please specifiy a noise level')
        return;
    end
    
end

disp('DMD error model is done')

%% DMD Learn Physics

clearvars -except param DMdata noise
dt = param.dt;
tlength = param.tlength;
tspan = param.tspan;
cutoff = param.cutoff;
order = param.order;
PD = param.PD;
g = param.g;
system = param.system;
r = param.r;
p = param.p;

noise = 0;
lowpass_filter = 0;

tStart = tic;
eval('DMD_discrep_IDphysics')
tEnd = toc(tStart);
    
DMdata.IDphysics.DMD.nonoise.xA = xA;
DMdata.IDphysics.DMD.nonoise.xB = xB;
DMdata.IDphysics.DMD.nonoise.xC = xC;
DMdata.IDphysics.DMD.nonoise.x = trainSet.X;
DMdata.IDphysics.DMD.nonoise.y = trainSet.Y;
DMdata.IDphysics.DMD.nonoise.xR = xR;
DMdata.IDphysics.DMD.nonoise.time = tEnd;

for noise = [0.001 0.01 0.1]
    
    clearvars -except param DMdata noise
    dt = param.dt;
    tlength = param.tlength;
    tspan = param.tspan;
    cutoff = param.cutoff;
    order = param.order;
    PD = param.PD;
    g = param.g;
    system = param.system;
    r = param.r;
    p = param.p;
        
    lowpass_filter = 1;
    
    tStart = tic;
    eval('DMD_discrep_IDphysics')
    tEnd = toc(tStart);
    
    if noise == 0.001
        DMdata.IDphysics.DMD.lownoise.xA = xA;
        DMdata.IDphysics.DMD.lownoise.xB = xB;
        DMdata.IDphysics.DMD.lownoise.xC = xC;
        DMdata.IDphysics.DMD.lownoise.x = trainSet.X;
        DMdata.IDphysics.DMD.lownoise.y = trainSet.Y;
        DMdata.IDphysics.DMD.lownoise.xR = xR;
        DMdata.IDphysics.DMD.lownoise.time = tEnd;
    elseif noise == 0.01
        DMdata.IDphysics.DMD.mednoise.xA = xA;
        DMdata.IDphysics.DMD.mednoise.xB = xB;
        DMdata.IDphysics.DMD.mednoise.xC = xC;
        DMdata.IDphysics.DMD.mednoise.x = trainSet.X;
        DMdata.IDphysics.DMD.mednoise.y = trainSet.Y;
        DMdata.IDphysics.DMD.mednoise.xR = xR;
        DMdata.IDphysics.DMD.mednoise.time = tEnd;
    elseif noise == 0.1
        DMdata.IDphysics.DMD.highnoise.xA = xA;
        DMdata.IDphysics.DMD.highnoise.xB = xB;
        DMdata.IDphysics.DMD.highnoise.xC = xC;
        DMdata.IDphysics.DMD.highnoise.x = trainSet.X;
        DMdata.IDphysics.DMD.highnoise.y = trainSet.Y;
        DMdata.IDphysics.DMD.highnoise.xR = xR;
        DMdata.IDphysics.DMD.highnoise.time = tEnd;
    else
        disp('Please specifiy a noise level')
        return;
    end
    
end

disp('DMD learn physics is done')

%% GPR Model Error

clearvars -except param DMdata noise
dt = param.dt;
tlength = param.tlength;
tspan = param.tspan;
cutoff = param.cutoff;
order = param.order;
PD = param.PD;
g = param.g;
system = param.system;
r = param.r;
p = param.p;

noise = 0;
lowpass_filter = 0;

tStart = tic;
eval('GPR_discrep_errorModel')
tEnd = toc(tStart);

DMdata.errorModel.GPR.nonoise.xA = xA;
DMdata.errorModel.GPR.nonoise.xB = xB;
DMdata.errorModel.GPR.nonoise.xC = xC;
DMdata.errorModel.GPR.nonoise.x = trainSet.X;
DMdata.errorModel.GPR.nonoise.y = trainSet.Y;
DMdata.errorModel.GPR.nonoise.xR = xR;
DMdata.errorModel.GPR.nonoise.time = tEnd;

for noise = [0.001 0.01 0.1]
    
    clearvars -except param DMdata noise
    dt = param.dt;
    tlength = param.tlength;
    tspan = param.tspan;
    cutoff = param.cutoff;
    order = param.order;
    PD = param.PD;
    g = param.g;
    system = param.system;
    r = param.r;
    p = param.p;
        
    lowpass_filter = 1;
    
    tStart = tic;
    eval('GPR_discrep_errorModel')
    tEnd = toc(tStart);
    
    if noise == 0.001
        DMdata.errorModel.GPR.lownoise.xA = xA;
        DMdata.errorModel.GPR.lownoise.xB = xB;
        DMdata.errorModel.GPR.lownoise.xC = xC;
        DMdata.errorModel.GPR.lownoise.x = trainSet.X;
        DMdata.errorModel.GPR.lownoise.y = trainSet.Y;
        DMdata.errorModel.GPR.lownoise.xR = xR;
        DMdata.errorModel.GPR.lownoise.time = tEnd;
    elseif noise == 0.01
        DMdata.errorModel.GPR.mednoise.xA = xA;
        DMdata.errorModel.GPR.mednoise.xB = xB;
        DMdata.errorModel.GPR.mednoise.xC = xC;
        DMdata.errorModel.GPR.mednoise.x = trainSet.X;
        DMdata.errorModel.GPR.mednoise.y = trainSet.Y;
        DMdata.errorModel.GPR.mednoise.xR = xR;
        DMdata.errorModel.GPR.mednoise.time = tEnd;
    elseif noise == 0.1
        DMdata.errorModel.GPR.highnoise.xA = xA;
        DMdata.errorModel.GPR.highnoise.xB = xB;
        DMdata.errorModel.GPR.highnoise.xC = xC;
        DMdata.errorModel.GPR.highnoise.x = trainSet.X;
        DMdata.errorModel.GPR.highnoise.y = trainSet.Y;
        DMdata.errorModel.GPR.highnoise.xR = xR;
        DMdata.errorModel.GPR.highnoise.time = tEnd;
    else
        disp('Please specifiy a noise level')
        return;
    end
    
end

disp('GPR error model is done')

%% GPR Learn Physics

clearvars -except param DMdata noise
dt = param.dt;
tlength = param.tlength;
tspan = param.tspan;
cutoff = param.cutoff;
order = param.order;
PD = param.PD;
g = param.g;
system = param.system;
r = param.r;
p = param.p;

noise = 0;
lowpass_filter = 0;

tStart = tic;
eval('GPR_discrep_IDphysics')
tEnd = toc(tStart);

DMdata.IDphysics.GPR.nonoise.xA = xA;
DMdata.IDphysics.GPR.nonoise.xB = xB;
DMdata.IDphysics.GPR.nonoise.xC = xC;
DMdata.IDphysics.GPR.nonoise.x = trainSet.X;
DMdata.IDphysics.GPR.nonoise.y = trainSet.Y;
DMdata.IDphysics.GPR.nonoise.xR = xR;
DMdata.IDphysics.GPR.nonoise.time = tEnd;

for noise = [0.001 0.01 0.1]
    
    clearvars -except param DMdata noise
    dt = param.dt;
    tlength = param.tlength;
    tspan = param.tspan;
    cutoff = param.cutoff;
    order = param.order;
    PD = param.PD;
    g = param.g;
    system = param.system;
    r = param.r;
    p = param.p;
        
    lowpass_filter = 1;
    
    tStart = tic;
    eval('GPR_discrep_IDphysics')
    tEnd = toc(tStart);
    
    if noise == 0.001
        DMdata.IDphysics.GPR.lownoise.xA = xA;
        DMdata.IDphysics.GPR.lownoise.xB = xB;
        DMdata.IDphysics.GPR.lownoise.xC = xC;
        DMdata.IDphysics.GPR.lownoise.x = trainSet.X;
        DMdata.IDphysics.GPR.lownoise.y = trainSet.Y;
        DMdata.IDphysics.GPR.lownoise.xR = xR;
        DMdata.IDphysics.GPR.lownoise.time = tEnd;
    elseif noise == 0.01
        DMdata.IDphysics.GPR.mednoise.xA = xA;
        DMdata.IDphysics.GPR.mednoise.xB = xB;
        DMdata.IDphysics.GPR.mednoise.xC = xC;
        DMdata.IDphysics.GPR.mednoise.x = trainSet.X;
        DMdata.IDphysics.GPR.mednoise.y = trainSet.Y;
        DMdata.IDphysics.GPR.mednoise.xR = xR;
        DMdata.IDphysics.GPR.mednoise.time = tEnd;
    elseif noise == 0.1
        DMdata.IDphysics.GPR.highnoise.xA = xA;
        DMdata.IDphysics.GPR.highnoise.xB = xB;
        DMdata.IDphysics.GPR.highnoise.xC = xC;
        DMdata.IDphysics.GPR.highnoise.x = trainSet.X;
        DMdata.IDphysics.GPR.highnoise.y = trainSet.Y;
        DMdata.IDphysics.GPR.highnoise.xR = xR;
        DMdata.IDphysics.GPR.highnoise.time = tEnd;
    else
        disp('Please specifiy a noise level');
        return;
    end
    
end

disp('GPR learn physics is done')

%% NN Model Error

clearvars -except param DMdata noise
dt = param.dt;
tlength = param.tlength;
tspan = param.tspan;
cutoff = param.cutoff;
order = param.order;
PD = param.PD;
g = param.g;
system = param.system;
r = param.r;
p = param.p;

noise = 0;
lowpass_filter = 0;

tStart = tic;
eval('NN_discrep_errorModel')
nnet.guis.closeAllViews()
tEnd = toc(tStart);

DMdata.errorModel.NN.nonoise.xA = xA;
DMdata.errorModel.NN.nonoise.xB = xB;
DMdata.errorModel.NN.nonoise.xC = xC;
DMdata.errorModel.NN.nonoise.x = trainSet.X;
DMdata.errorModel.NN.nonoise.y = trainSet.Y;
DMdata.errorModel.NN.nonoise.xR = xR;
DMdata.errorModel.NN.nonoise.time = tEnd;

for noise = [0.001 0.01 0.1]
    
    clearvars -except param DMdata noise
    dt = param.dt;
    tlength = param.tlength;
    tspan = param.tspan;
    cutoff = param.cutoff;
    order = param.order;
    PD = param.PD;
    g = param.g;
    system = param.system;
    r = param.r;
    p = param.p;
        
    lowpass_filter = 1;
    
    tStart = tic;
    eval('NN_discrep_errorModel')
    tEnd = toc(tStart);
    
    if noise == 0.001
        DMdata.errorModel.NN.lownoise.xA = xA;
        DMdata.errorModel.NN.lownoise.xB = xB;
        DMdata.errorModel.NN.lownoise.xC = xC;
        DMdata.errorModel.NN.lownoise.x = trainSet.X;
        DMdata.errorModel.NN.lownoise.y = trainSet.Y;
        DMdata.errorModel.NN.lownoise.xR = xR;
        DMdata.errorModel.NN.lownoise.time = tEnd;
    elseif noise == 0.01
        DMdata.errorModel.NN.mednoise.xA = xA;
        DMdata.errorModel.NN.mednoise.xB = xB;
        DMdata.errorModel.NN.mednoise.xC = xC;
        DMdata.errorModel.NN.mednoise.x = trainSet.X;
        DMdata.errorModel.NN.mednoise.y = trainSet.Y;
        DMdata.errorModel.NN.mednoise.xR = xR;
        DMdata.errorModel.NN.mednoise.time = tEnd;
    elseif noise == 0.1
        DMdata.errorModel.NN.highnoise.xA = xA;
        DMdata.errorModel.NN.highnoise.xB = xB;
        DMdata.errorModel.NN.highnoise.xC = xC;
        DMdata.errorModel.NN.highnoise.x = trainSet.X;
        DMdata.errorModel.NN.highnoise.y = trainSet.Y;
        DMdata.errorModel.NN.highnoise.xR = xR;
        DMdata.errorModel.NN.highnoise.time = tEnd;
    else
        disp('Please specifiy a noise level');
        return;
    end
    
end

disp('NN error model is done')


%% NN Learn Physics

clearvars -except param DMdata noise
dt = param.dt;
tlength = param.tlength;
tspan = param.tspan;
cutoff = param.cutoff;
order = param.order;
PD = param.PD;
g = param.g;
system = param.system;
r = param.r;
p = param.p;

noise = 0;
lowpass_filter = 0;

tStart = tic;
eval('NN_discrep_IDphysics')
nnet.guis.closeAllViews()
tEnd = toc(tStart);

DMdata.IDphysics.NN.nonoise.xA = xA;
DMdata.IDphysics.NN.nonoise.xB = xB;
DMdata.IDphysics.NN.nonoise.xC = xC;
DMdata.IDphysics.NN.nonoise.x = trainSet.X;
DMdata.IDphysics.NN.nonoise.y = trainSet.Y;
DMdata.IDphysics.NN.nonoise.xR = xR;
DMdata.IDphysics.NN.nonoise.time = tEnd;

for noise = [0.001 0.01 0.1]
    
    clearvars -except param DMdata noise
    dt = param.dt;
    tlength = param.tlength;
    tspan = param.tspan;
    cutoff = param.cutoff;
    order = param.order;
    PD = param.PD;
    g = param.g;
    system = param.system;
    r = param.r;
    p = param.p;
        
    lowpass_filter = 1;
    
    tStart = tic;
    eval('NN_discrep_IDphysics')
    tEnd = toc(tStart);
    
    if noise == 0.001
        DMdata.IDphysics.NN.lownoise.xA = xA;
        DMdata.IDphysics.NN.lownoise.xB = xB;
        DMdata.IDphysics.NN.lownoise.xC = xC;
        DMdata.IDphysics.NN.lownoise.x = trainSet.X;
        DMdata.IDphysics.NN.lownoise.y = trainSet.Y;
        DMdata.IDphysics.NN.lownoise.xR = xR;
        DMdata.IDphysics.NN.lownoise.time = tEnd;
    elseif noise == 0.01
        DMdata.IDphysics.NN.mednoise.xA = xA;
        DMdata.IDphysics.NN.mednoise.xB = xB;
        DMdata.IDphysics.NN.mednoise.xC = xC;
        DMdata.IDphysics.NN.mednoise.x = trainSet.X;
        DMdata.IDphysics.NN.mednoise.y = trainSet.Y;
        DMdata.IDphysics.NN.mednoise.xR = xR;
        DMdata.IDphysics.NN.mednoise.time = tEnd;
    elseif noise == 0.1
        DMdata.IDphysics.NN.highnoise.xA = xA;
        DMdata.IDphysics.NN.highnoise.xB = xB;
        DMdata.IDphysics.NN.highnoise.xC = xC;
        DMdata.IDphysics.NN.highnoise.x = trainSet.X;
        DMdata.IDphysics.NN.highnoise.y = trainSet.Y;
        DMdata.IDphysics.NN.highnoise.xR = xR;
        DMdata.IDphysics.NN.highnoise.time = tEnd;
    else
        disp('Please specifiy a noise level')
        return;
    end
    
end

disp('NN learn physics is done')

%% Save results

save([system,'_All_Results.mat'],'DMdata')
