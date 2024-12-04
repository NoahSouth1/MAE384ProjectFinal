%% Section 3
% Q1 - True Data Set (30 day model)
for j = 1:2 
    %j = 1 = Q1 - True Data Set (30 day model)
    %j = 2 = Q2 - Estimate via Linear Least Sq (30 day model)
    %j = 3 = Q3 - Estimate via Linear Least Sq (10 day model)

    if j == 1
        h = 1; 
        S0 = 990; 
        I0 = 10; 
        R0 = 0; 
        T = 30;
        time_vec = 0:h:T; 

        steps = length(time_vec); 
        
        % Parameter sets
        beta_true = 0.3; 
        gam_true = 0.1;
        
        figure(j);
        
        % True data set
        [S, I, R] = runge_kutta(S0, I0, R0, beta_true, gam_true, h, steps);
        plot(time_vec, S, 'b', time_vec, I, 'r', time_vec, R, 'g');
        title('True Data Set I(t)');
        xlabel('Time');
        ylabel('Population');
        legend('Susceptible', 'Infected', 'Recovered');
        grid on;
    
    else
        h = 1; 
        S0 = 990; 
        I0 = 10; 
        R0 = 0; 
        T = 10;
        time_vec = 1:h:T; 
        
        steps = length(time_vec); 
        
        % Parameter sets
        beta_true = 0.3; 
        gam_true = 0.1;
        N = 1000;

        k = ((beta_true*S0)/N) - gam_true;

        figure(j);
        
        %for z = 1:steps
        %    I0_est(z) = esp(log(I(z)) - (k*time_vec(1,z)))
        %end

        x_input = time_vec.';
        y_output = I(2:11);
        x2 = x_input.^2;
        xy = x_input .* y_output;
        x_avg = mean(x_input);
        y_avg = mean(y_output);

        a1 = (N*sum(xy)) - (sum(x_input))*(sum(y_output));
        a1 = a1 / (N*sum(x2) - (sum(x_input)^2));

        a0 = y_avg - (a1*x_avg);

        y_graph = a0 + (time_vec.*a1);
        plot(x_input, y_graph, 'r');
        hold on
        scatter(x_input,y_output)

        I0_est = a0

        

        
        
    
    end
end


%% end of code, begin of function library

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
