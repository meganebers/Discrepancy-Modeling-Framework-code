function [RMSE_plato, RMSE_aug, perChange] = Compare_RSME(xA,xB,xC,n,time)
%butterfilt computes Butterworth filter parameters and filters are the
%precribed frequency
% Required INPUTS first
% xA - true data
% xB - augmented data
% xC - approximate data
% n - number of many states (columns)
% 
% Optional INPUTS
% time - Time column; only induce if Data does not include time
%
% OUTPUT
% RMSE_plato - root mean squared error of true data vs approximate data
% RMSE_aug - root mean squared error of true data vs augmented data (e.g., discrepancy + approximate model
% perChange - percent change in RMSE by augmenting the approximate model
%
% Megan Ebers
% University of Washington - Autumn 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check for good data
[rd,cd]=size(xA); % Get data dimensions
if rd<cd % Assume that the short dimension denotes states. Make this the column dimension
    xA=xA';
end

[rd,cd]=size(xB); % Get data dimensions
if rd<cd % Assume that the short dimension denotes states. Make this the column dimension
    xB=xB';
end

[rd,cd]=size(xC); % Get data dimensions
if rd<cd % Assume that the short dimension denotes states. Make this the column dimension
    xC=xC';
end

if ~exist('time') % Assume that time is in the first column of 'data'. create 'time' variabel
    time=xA(:,1);
    tt = length(time);
else % Ensure time is a column vector
    tt = time;
    %[rt,ct]=size(time);
    %if rt<ct
    %    time=time';
    %end
end

if n>1
    for i = 1:n
        RMSE_plato(:,i) = sqrt(mean((xA(1:tt,i) - xC(1:tt,i)).^2));
        RMSE_aug(:,i) = sqrt(mean((xA(1:tt,i) - xB(1:tt,i)).^2));
        perChange(i) = (RMSE_plato(i)-RMSE_aug(i))/(RMSE_plato(i))*100;
    end
else
    RMSE_plato = sqrt(mean((xA(1:tt,1) - xC(1:tt,1)).^2));
    RMSE_aug = sqrt(mean((xA(1:tt,1) - xB(1:tt,1)).^2));
    perChange = (RMSE_plato-RMSE_aug)/(RMSE_plato)*100;
end

disp(['RMSE of Platonic state space error: ',num2str(RMSE_plato)])
disp(['RMSE of Augemented state space error: ',num2str(RMSE_aug)])
disp(['Percent decrease in RMSE: ',num2str(perChange),'%'])
