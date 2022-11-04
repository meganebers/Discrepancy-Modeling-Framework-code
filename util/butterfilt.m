function [ dataFilt ] = butterfilt( data,cutoff,order, pass,time,plt  )
%butterfilt computes Butterworth filter parameters and filters are the
%precribed frequency
% Required INPUTS first
% data - Data to filter
% cutoff - cutoff frequecy (Hz)
% order - nth order Butterworth filter
% pass - 'High','Low' pass filter
% plt - plotting bolean ('on','off')
% 
% Optional INPUTS
% time - Time column; only inluce if Data does not include time
%
% OUTPUT
% datafilt - filtered data
%
% Michael Rosenberg
% University of Washington - Summer 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check for good data
[rd,cd]=size(data); % Get data dimensions
if rd<cd % Assume that the short dimension denotes states. Make this the column dimension
    data=data';
end

if ~exist('time') % Assume that time is in the first column of 'data'. create 'time' variabel
    time=data(:,1);
else % Ensure time is a column vector
    [rt,ct]=size(time);
    if rt<ct
        time=time';
    end
end

if ~all(diff(time)>0) % Ensure monotonic time
    error('Variable "time" must be monotonically-increasing. Ensure that the first row or column of data represents time and ensure that you are inputting a [time x state] array.');
end

% Check that other inputs are OK
if mod(order,floor(order))~=0;   error('Butterworth filter order must be a whole number.'); end
if	~strcmp(pass,'High') && ~strcmp(pass,'Low');   error('"pass" must be either the string "High" or "Low".'); end
%% Compute Butterworth filter inputs
nyq=1/mean(diff(time))/2; % Nyquist frequency
k=cutoff/nyq; % Filter "gain"

%% Filter
[b,a]=butter(order,k,pass); % Butterworth filter parameters
dataFilt=filtfilt(b,a,data); % Filter data 
% Plot so user can verify correct signal
if strcmp(plt,'on')
figure
    hold on
    col=randi([1,length(data(1,:))]);
    plot(data(:,col),'k'); plot(dataFilt(:,col),'b');
    legend('Original signal','Filtered signal')
    xlabel('Time'); ylabel(strcat('Data column #',num2str(col)));
end
end