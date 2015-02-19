function [gaps,temperatures] = main_superconductor_critical()
    pos_len = 8;
    erg_len = 8;
    cutoff = 20;
    
    % Create a superconductor with N₀λ = 0.2, T=0
    s = Superconductor( linspace(0,1,pos_len), [linspace(0,1.5,erg_len) linspace(1.6,cutoff,erg_len)], 300, 1, 0.2 );
    s.scaling = 0.2;
    s.temperature = 1e-16;
    
    % Bootstrap the state by solving the equations self-consistently at T=0
    for n=1:3
        fprintf(':: Initializing the superconductor at absolute zero... [ iteration %.f ]\n', n);
        s.update;
    end
    
    % Save the T=0 results to the output vector
    temperatures = [0];
    gaps         = [mean(s.gap)];

    % Keep increasing the temperature and measuring the gap until critical
    while ~s.critical
        % Increase the temperature of the system, and update the state
        s.temperature = s.temperature + 0.1;
        s.update;
        
        % Save the current temperature and gap to output vectors
        temperatures(end+1) = s.temperature;
        gaps(end+1) = mean(s.gap);
    end
    
    % Plot the results
    figure;
    title('Plot of temperature vs. superconducting gap');
    plot(temperatures, gaps);
    xlabel('Temperature');
    ylabel('Superconducting gap');
end