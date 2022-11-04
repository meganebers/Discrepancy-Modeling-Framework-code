function plotGPR(xtest,ytest,ypred,yci)
figure,
plot(xtest,ytest,'r.');
hold on;
plot(xtest,ypred);
hold on;
plot(xtest,yci(:,1),'k:');
hold on;
plot(xtest,yci(:,2),'k:');
xlabel('x');
ylabel('y');
legend('95% CI','Measurements','GPR predictions')
hold off




