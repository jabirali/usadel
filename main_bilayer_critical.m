function main_bilayer_critical()
    % Create a superconductor and ferromagnet
    s = Superconductor( linspace(0,1,32), [linspace(0.1,1.5,0.1) linspace(5,40,5)], 300, 1, 0.2);
    f = Ferromagnet(    linspace(0,1,32), [linspace(0.1,1.5,0.1) linspace(5,40,5)], 300, 1);

    % Make the interface between the materials transparent
    s.interface_right = 1;
    f.interface_left  = 1;

    % Place to store the output
    temperatures = [];
    gaps         = [];
    
    while ~s.critical
        % Exchange boundary conditions
        s.boundary_right(:) = f.states(1,  :);
        f.boundary_left(:)  = s.states(end,:);
        
        % Update the system
        s.update;
        f.update;
        
        % Save the current gap and temperature
        temperatures(end+1) = s.temperature;
        gaps(end+1)         = mean(abs(s.gap));
        
        % Increase the temperature of the system
        s.temperature = s.temperature + 0.1;
        fprintf(':: INCREASING TEMPERATURE TO %f.\n', s.temperature);
    end
end