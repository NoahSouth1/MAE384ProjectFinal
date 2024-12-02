for i=1:2
    h = 0.1; 
    S0 = 990; 
    I0 = 10; 
    R0 = 0; 
    T = 30; 
    time_vec = 0:h:T;
    steps = length(time_vec);
    beta = 0.3;
    beta_var = zeros(1,steps);
    if i==1
        for k = 1 : steps
            beta_var(k) = beta * (1 + 5 * sin(2 * pi * (365/365) * time_vec(k)));
        end 
    else
        for k = 1 : steps
            beta_var(k) = beta * (1 + 5 * sin(2 * pi * (100/365) * time_vec(k)));
        end 
    end
    gamma = 0.1;
    
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
        dS1 = -beta_var(t) * S(t) * I(t) / N;
        dI1 = beta_var(t) * S(t) * I(t) / N - gamma * I(t);
        dR1 = gamma * I(t);
    
        % k2
        dS2 = -beta_var(t) * (S(t) + h*dS1/2) * (I(t) + h*dI1/2) / N;
        dI2 = beta_var(t) * (S(t) + h*dS1/2) * (I(t) + h*dI1/2) / N - gamma * (I(t) + h*dI1/2);
        dR2 = gamma * (I(t) + h*dI1/2);
    
        % k3
        dS3 = -beta_var(t) * (S(t) + h*dS2/2) * (I(t) + h*dI2/2) / N;
        dI3 = beta_var(t) * (S(t) + h*dS2/2) * (I(t) + h*dI2/2) / N - gamma * (I(t) + h*dI2/2);
        dR3 = gamma * (I(t) + h*dI2/2);
    
        % k4
        dS4 = -beta_var(t) * (S(t) + h*dS3) * (I(t) + h*dI3) / N;
        dI4 = beta_var(t) * (S(t) + h*dS3) * (I(t) + h*dI3) / N - gamma * (I(t) + h*dI3);
        dR4 = gamma * (I(t) + h*dI3);
    
        % Update SIR
        S(t+1) = S(t) + h * (dS1 + 2*dS2 + 2*dS3 + dS4) / 6;
        I(t+1) = I(t) + h * (dI1 + 2*dI2 + 2*dI3 + dI4) / 6;
        R(t+1) = R(t) + h * (dR1 + 2*dR2 + 2*dR3 + dR4) / 6;
    end
    
    %Fast Fourier Transform
    Sf=fft(S);
    If=fft(I);
    Rf=fft(R);
    %Spectrum
    Is=zeros(1,(steps-1)/2);
    k=1;
    while k<=(steps+1)/2
        Is(k)=abs(If(k));
        k=k+1;
    end
    f=(0:(steps-1)/2)/T;
    %Plot Results
    if i==1
    figure(1)
        subplot(3,1,1); %SIR Model
            plot(time_vec,S,'b', time_vec,I,'r', time_vec,R,'g');
            title('SIR Model');
            xlabel('Time');
            ylabel('Population');
            legend('Susceptible', 'Infected', 'Recovered');
        subplot(3,1,2); %Transmission Rate
            plot(time_vec,beta_var);
            title('Transmission Rate');
            xlabel('time');
            ylabel('Magnitude');
        subplot(3,1,3) %Spectrum
            plot(f,Is);
            title('Discrete Fourier Transform');
            xlabel('Frequency');
            ylabel('|Fourier coefficients|');
    else
        figure(2)
        subplot(3,1,1); %SIR Model
            plot(time_vec,S,'b', time_vec,I,'r', time_vec,R,'g');
            title('SIR Model');
            xlabel('Time');
            ylabel('Population');
            legend('Susceptible', 'Infected', 'Recovered');
        subplot(3,1,2); %Transmission Rate
            plot(time_vec,beta_var);
            title('Transmission Rate');
            xlabel('time');
            ylabel('Magnitude');
        subplot(3,1,3) %Spectrum
            plot(f,Is);
            title('Discrete Fourier Transform');
            xlabel('Frequency');
            ylabel('|Fourier coefficients|');
    end
end