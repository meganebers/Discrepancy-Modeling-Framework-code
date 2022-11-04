function [ P, T ] = estimate_phaseBivariate( X1,X2,dt,plt, filterSpec, X1_mean, X2_mean, X1_sd, X2_sd)
%estimate_phaseBivariate generates phase estimates as the phase and angle defined by two input states
%
% - The phase corresponds to the angle of the portrait X1+i*X2
% - To wrap P, use P = mod(P, 2*pi);
% - For mathematical convenience or if investigating changes in phase, make phase monotonically-increasing. I
% don't think it's super beneficial or necessary for all applications.
%
% Steps
% 0. Define some convenient parameters
% 1. Scale variables to z-scores
% 2. Compute phase and filter to ensure monotonicity
% 3. Get period
% 4. Plots
%
% INPUTS
% X1 - array [N x 1] Time history (N samples) of the first phase variable
% X2 - array [N x 1] (optional) Time history (N samples) of the first
% second variable. If X2 is empty, X2 is defined as the derivative of X1 with time step dt
% dt - scalar Time step
% plt - string to turn plotting 'on',or 'off' 
% filterSpec - 'on' or 'off' specify whether to low-pass filter or not to achieve monotonically-increasing phase
% Use The following if you want to ensure consistent phase estimates across subsets of the gait cycle
%   X1_mean - scalar mean of phase variable X1; 
%   X2_mean - scalar mean of phase variable X2; 
%   X1_sd - scalar standard deviation of phase variable X1;
%   X2_sd- scalar standard deviation of phase variable X2;
%
% OUTPUTS
% P - [N x 1] Unwrapped phase estimate for the input dataset
% T - [N x 1] Cycle period
%
% Michael Rosenberg - University of Washington - 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Preliminary
[N1,~] = size(X1);
if isempty(X2)
    X2 = gradient(X1,dt);
end
[N2,~] = size(X2);
if N1 ~= N2
    error('Arrays must be the same size!'); 
end

% Ste defaults
if ~exist('filterSpec','var') % To filter or not to filter...
    filterSpec = 'on';
end
%% 1. Studentize  - De-mean & scale to unit variance
% Store originals
X1_old=X1;
X2_old=X2;
% Studentize
if ~exist('X1_mean','var')
    X1=(X1-mean(X1))./std(X1);
    X2=(X2-mean(X2))./std(X2);
else % Must de-mean with provided values
    X1=(X1-X1_mean) ./ X1_sd;
    X2=(X2-X2_mean) ./ X2_sd;
end

%% 2. Compute phase and unwrap
Pwrapped = mod(angle(complex(X1,X2)), 2*pi); % Wrapped phase [0,2*pi)
P=unwrap(Pwrapped); % this is out phase estimate (in radians) - need to unwrap to filter

% Check for monotonicity & increasing value
if sum(P) < 0 % Ensure that phase increases in time
    P=-P;
end

% If desired, low-pass filter the data until phase is monotonically-increasing
dp=diff(P);
if ~isempty(find(dp<0)) && ~exist('X1_mean','var') && strcmp(filterSpec, 'on')
% %     disp('WARNING: Phase is not monotonically increasing!');
    % Attempt to smooth the phase using a low-pass filter
% %     disp('Attempting to smooth phase using a low-pass filter')
    LP = 6;
    while ~isempty(find(dp<0)) && LP>.25 % iteratively reduce LP until monotonicity is achieved
        LP = LP*.75;
        P_monotonic = butterfilt(P,LP,4,'Low',linspace(dt,N1*dt,N1),'off'); 
        dp = diff(P_monotonic);
    end
    P_monotonic=butterfilt(P,LP,4,'Low',linspace(dt,N1*dt,N1),'off');
% %     disp(strcat('Phase is monotonically increasing once filtered at LP=',num2str(LP),'Hz.'))
    Pwrapped=mod(P_monotonic,2*pi)-pi;
else % If we don't care about phase increasing monotonically... 
    P = P; % Just so you know nothing's happening
end
%% 3. Get period - rarely necessary
ct = 1;
xition = [];
for n = 1:N1-1
    if Pwrapped(n) > Pwrapped(n+1)% This is all we need, since phase is monotonically increasing (not elegant)
        xition(ct) = n;
        ct = ct+1;
    end
end
if length(xition>1)
    T = mean(diff(xition)); % Period
    Ts = std(diff(xition)); % Period SD
% %     disp(strcat('The average +/- SD period of the phase is:...',num2str(T),'+/-',num2str(Ts),'...Frames'));
else
    T = [];Ts = [];
end
%% 4. Plots 
% Phase portrait
if strcmp(plt,'on')
    figure
    plot(X1,X2,'k.')
    xlabel('X1');ylabel('X2')

% Phase
    figure
    plot(Pwrapped,'k')
    xlabel('Data point');ylabel('Phase (rad)')
    
% Inputs in phase domain
    figure
    subplot(1,2,1)
    plot(Pwrapped,X1_old,'k.'); 
    xlabel('Phase (rad)');ylabel('X1')
    subplot(1,2,2)
    plot(Pwrapped,X2_old,'k.')
    xlabel('Phase (rad)');ylabel('X2')
end

end



