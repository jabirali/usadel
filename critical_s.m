function [gaps,temperatures] = critical_s()
% TODO: Switch to scaling by correlation length?
    pos_len = 16;
    erg_len = 8;
    cutoff = 20;
    
    temperatures = [];
    gaps         = [];
    parfor n=1:24    
        % Create a superconductor at a new temperature
        s = Superconductor( linspace(0,1,pos_len), [linspace(0,1.5,erg_len) linspace(1.6,cutoff,erg_len)] );
        s.scaling     = 0.0001;
        s.temperature = n;
        
        % Update the state of the superconductor
        s.gap_update;
        for i=1:3
            s.update;
        end
 
        % Save the current temperature and gap to output vectors
        temperatures(n) = s.temperature;
        gaps(n)         = mean(s.gap);

        % Can't break parallel for loop, so let's just return info instead
        if ~s.critical
           fprintf('Note: Critical at T=%f.', s.temperature);
        end        
    end
    
    % Plot the results
    figure;
    title('Plot of temperature vs. superconducting gap');
    plot(temperatures, gaps);
    xlabel('Temperature');
    ylabel('Superconducting gap');
end