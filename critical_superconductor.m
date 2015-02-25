% This script calculates the superconducting gap as a function of temperature.

% Define the simulation parameters
output       = 'output/critical_superconductor/';
positions    = linspace(0, 1, 64);
temperatures = linspace(0, 1, 20);
energies     = [linspace(1e-5, 1.5, 20) linspace(1.55,75,10)];
%energies     = [0.00001:0.01:0.90, 0.91:0.005:1.09, 1.1:0.01:1.5, 2:1:10, 15:10:75];
thouless     = 4;
strength     = 0.2;

% Create a superconductor based on the parameters above
s = Superconductor(positions, energies, thouless, strength);

% Create the output directory if it doesn't already exist
mkdir(output);

% Define variables used to store results
gap  = 0;                           % Current superconducting gap
gaps = zeros(size(temperatures));   % All the calculated gaps

for n=1:length(temperatures)
    % Progress information
    fprintf('\n:: PROGRAM:        [ %3d / %3d ]  T = %d\n', n, length(temperatures), temperatures(n));
    
    % Increase the superconductor temperature
    s.temperature = temperatures(n);
    
    % Keep updating the internal state of the superconductor until the gap
    % converges (i.e. less than 0.1% change between iterations)
    while abs(1 - s.gap(1/2)/gap) > 1e-3
        gap = s.gap(1/2);
        s.update;
    end
    
    % Store the current result in the vector 'gaps'
    gaps(n) = s.gap(1/2);
end

% Save the results to file
save([output, 'output.mat'], 'temperatures', 'gaps')
    
% Plot the results
figure;
title('Plot of temperature vs. superconducting gap');
plot(temperatures, gaps);
xlabel('Temperature');
ylabel('Superconducting gap');