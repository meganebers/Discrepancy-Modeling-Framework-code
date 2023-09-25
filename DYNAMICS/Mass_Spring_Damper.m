%% NL mass spring damper with forcing

T = 0:0.01:20;
S = [0.5 0];

C = [1+0.1 0]; % Modified C matrix with measurement bias

[t, y_L] = ode45(@odefun_lin, T, S); % linear model
y_L_bias = y_L * C'; % Apply measurement bias to the output
y_L_bias =y_L_bias + noise*randn(size(y_L_bias)); % Apply measurement bias to the output

[t, y_NL] = ode45(@odefun_NL, T, S); % nonlinear model
y_NL_bias = y_NL * C'; % Apply measurement bias to the output
y_NL_bias =y_NL_bias + noise*randn(size(y_L_bias)); % Apply measurement bias to the output

%{
figure, plot(t,y_L, '-', t,y_NL, '--', 'Linewidth', [2])
legend({'Position', 'Speed'});
ylabel('Position / Speed [m / m/s]')
xlabel('Time [s]')
title('Mass-Spring-Damper system (linear vs nonlinear)');

figure, plot(t,y_L, '-', t, y_L_bias, '--', 'Linewidth', [2])
legend({'Position', 'Speed'});
ylabel('Position / Speed [m / m/s]')
xlabel('Time [s]')
title('Mass-Spring-Damper system (linear vs measurement bias)');

figure, plot(t,y_NL, '-', t, y_NL_bias, '--', 'Linewidth', [2])
legend({'Position', 'Speed'});
ylabel('Position / Speed [m / m/s]')
xlabel('Time [s]')
title('Mass-Spring-Damper system (nonlinear vs measurement bias)');
%}

function dydt = odefun_lin(~,x)
    F = 1*x(1)^3 ;%- 1*x(1)^2 + 1*x(1);
    dydt = [x(2);  F - (0.5/2)*x(2) - (2/2)*x(1)];
end

function dydt = odefun_NL(~,x)
    F = 1*x(1)^3 ;%- 1*x(1)^2 + 1*x(1);
    dydt = [x(2); F - (0.5/2)*x(2) - (2/2)*(x(1) - 0.5*x(1).^3)];
end