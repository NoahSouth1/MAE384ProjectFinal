h = 1; 
S0 = 990; 
I0 = 10; 
R0 = 0; 
T = 100; 
time_vec = 0:h:T; 
steps = length(time_vec); 

% Parameter sets
beta_influenza = 0.3; 
gam_influenza = 0.1;

beta_covid = 1.0; 
gam_covid = 0.1;

beta_measles = 2.0; 
gam_measles = 0.2;

figure;

% Seasonal Influenza
[S, I, R] = runge_kutta(S0, I0, R0, beta_influenza, gam_influenza, h, steps);
subplot(3, 1, 1);
plot(time_vec, S, 'b', time_vec, I, 'r', time_vec, R, 'g');
title('Seasonal Influenza');
xlabel('Time');
ylabel('Population');
legend('Susceptible', 'Infected', 'Recovered');
grid on;

% COVID-19
[S, I, R] = runge_kutta(S0, I0, R0, beta_covid, gam_covid, h, steps);
subplot(3, 1, 2);
plot(time_vec, S, 'b', time_vec, I, 'r', time_vec, R, 'g');
title('COVID-19');
xlabel('Time');
ylabel('Population');
legend('Susceptible', 'Infected', 'Recovered');
grid on;

% Measles
[S, I, R] = runge_kutta(S0, I0, R0, beta_measles, gam_measles, h, steps);
subplot(3, 1, 3);
plot(time_vec, S, 'b', time_vec, I, 'r', time_vec, R, 'g');
title('Measles');
xlabel('Time');
ylabel('Population');
legend('Susceptible', 'Infected', 'Recovered');
grid on;

% Runge-Kutta
function [S, I, R] = runge_kutta(S0, I0, R0, beta, gamma, h, steps)
    S = zeros(steps, 1);
    I = zeros(steps, 1);
    R = zeros(steps, 1);

    % initial conditions
    S(1) = S0;
    I(1) = I0;
    R(1) = R0;

    % Total pop
    N = S0 + I0 + R0;

    for t = 1:steps-1
        % k1
        dS1 = -beta * S(t) * I(t) / N;
        dI1 = beta * S(t) * I(t) / N - gamma * I(t);
        dR1 = gamma * I(t);

        % k2
        dS2 = -beta * (S(t) + h*dS1/2) * (I(t) + h*dI1/2) / N;
        dI2 = beta * (S(t) + h*dS1/2) * (I(t) + h*dI1/2) / N - gamma * (I(t) + h*dI1/2);
        dR2 = gamma * (I(t) + h*dI1/2);

        % k3
        dS3 = -beta * (S(t) + h*dS2/2) * (I(t) + h*dI2/2) / N;
        dI3 = beta * (S(t) + h*dS2/2) * (I(t) + h*dI2/2) / N - gamma * (I(t) + h*dI2/2);
        dR3 = gamma * (I(t) + h*dI2/2);

        % k4
        dS4 = -beta * (S(t) + h*dS3) * (I(t) + h*dI3) / N;
        dI4 = beta * (S(t) + h*dS3) * (I(t) + h*dI3) / N - gamma * (I(t) + h*dI3);
        dR4 = gamma * (I(t) + h*dI3);

        % Update SIR
        S(t+1) = S(t) + h * (dS1 + 2*dS2 + 2*dS3 + dS4) / 6;
        I(t+1) = I(t) + h * (dI1 + 2*dI2 + 2*dI3 + dI4) / 6;
        R(t+1) = R(t) + h * (dR1 + 2*dR2 + 2*dR3 + dR4) / 6;
    end
end

