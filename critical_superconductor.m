% This script calculates the superconducting gap as a function of temperature.

% Define the simulation parameters
output       = 'output/critical_superconductor/';
positions    = linspace(0, 1, 100);
temperatures = linspace(0, 1, 20);
energies     = [linspace(1e-5, 1.5, 48) linspace(1.55,75,12)];
thouless     = 1/16;
strength     = 0.2;

% Create a superconductor based on the parameters above
s = Superconductor(positions, energies, thouless, strength);

% Create the output directory (if it doesn't already exist)
mkdir(output);

% Log the program output to a file
diary([output, 'log.txt']);

% Define variables used to store results
gap  = 1;                          % Current superconducting gap
gaps = zeros(size(temperatures));  % All the calculated gaps

for n=1:length(temperatures)    
    % Status information
    fprintf('\n:: PROGRAM:        [ %3d / %3d ]  T = %.6f, gap = %.6f\n\n',  n, length(temperatures), temperatures(n), gap);

    % Increase the superconductor temperature, and update the gap and state
    s.temperature = temperatures(n);
    s.update_gap;
    s.update;
    
    % Keep updating the internal state of the superconductor until the gap
    % converges (i.e. less than 1% change between iterations)
    while (abs(1 - s.mean_gap/gap) > 1e-2) && ~s.critical
        % Status information
        fprintf('\n:: PROGRAM:        [ %3d / %3d ]  T = %.6f, gap = %.6f\n                   Gap changed by %.2f%%. Recalculating...\n\n',  n, length(temperatures), temperatures(n), s.mean_gap, 100*abs(1-s.mean_gap/gap));
            
        % Update the superconductor state
        gap = s.mean_gap;
        s.update;
    end
    
    % Store the current result in the vector 'gaps'
    gaps(n) = s.mean_gap;
end

% Save the results to file
save([output, 'output.mat'], 'temperatures', 'gaps');

% Disable the logging
diary off;

% Plot the results
figure;
title('Plot of temperature vs. superconducting gap');
plot(temperatures, gaps);
xlabel('Temperature');
ylabel('Superconducting gap');