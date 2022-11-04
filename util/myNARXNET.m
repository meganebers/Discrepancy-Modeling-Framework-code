function a = myNARXNET(in, out)

X = tonndata(in,false,false);
T = tonndata(out,false,false);

trainFcn = 'trainlm';  % Levenberg-Marquardt backpropagation.

inputDelays = 1:2;
feedbackDelays = 1:2;
hiddenLayerSize = [10,10,10];

net = narxnet(inputDelays,feedbackDelays,hiddenLayerSize,'open',trainFcn);

net.layers{1}.transferFcn = 'radbasn';
net.layers{2}.transferFcn = 'tansig';
net.layers{3}.transferFcn = 'tansig';

net.inputs{1}.processFcns = {'mapminmax','mapstd'};
net.inputs{2}.processFcns = {'mapminmax','mapstd'};

[x,xi,ai,t] = preparets(net,X,{},T);

net.trainParam.max_fail = 10;
net.trainParam.goal = 1e-10;

net.performParam.regularization = 0.5;

net.trainParam.epochs = 1000;
net.trainParam.min_grad = 1e-10;

net.divideFcn = 'dividerand';  % Divide data randomly
net.divideMode = 'time';  % Divide up every sample

net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

net.performFcn = 'mse';  % Mean Squared Error

net.plotFcns = {'plotperform','plottrainstate', 'ploterrhist', ...
    'plotregression', 'plotresponse', 'ploterrcorr', 'plotinerrcorr'};

% Train the Network
[net,tr] = train(net,x,t,xi,ai);
y = net(x,xi,ai);
e = gsubtract(t,y);
performance = perform(net,t,y);

%Switch the network to closed loop form and retrain
net = closeloop(net);
[x,xi,ai,t] = preparets(net,X,{},T);
[net,t] = train(net,x,t,xi,ai);
y = net(x,xi,ai);
a = net;    

end