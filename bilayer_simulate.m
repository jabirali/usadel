function bilayer_simulate()
    % Create a superconductor and ferromagnet
    s = Superconductor( linspace(0,1,32), [linspace(0.1,1.5,0.1) linspace(5,40,5)], 300, 0.2, 1);
    f = Ferromagnet(    linspace(0,1,32), [linspace(0.1,1.5,0.1) linspace(5,40,5)], 300, 1, 0, 1, pi/4);

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
        
        % Update the system using two parallel workers
        spmd
            switch labindex
                case 1
                    s = update_superconductor(s);
                case 2
                    f = update_ferromagnet(f):
            end
        end
        
        % TODO: Plot current singlet/triplet distribution in entire system
       
        % Save the current temperature and gap to the output vectors
        temperatures(end+1) = s.temperature;
        gaps(end+1)         = mean(abs(s.gap));
        
        % Increase the temperature of the system
        s.temperature = s.temperature + 0.1;
    end
end

function s = update_superconductor(input)
    s = input;
    s.update;
end

function f = update_ferromagnet(input)
    f = input;
    f.update;
end