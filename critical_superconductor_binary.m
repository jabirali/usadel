% This script calculates the critical temperature of a bulk superconductor,
% by performing a binary search for the temperature where the gap vanishes
% numerically. The result should be numerically one in the given unit
% system, so this script can be used to calibrate simulation parameters.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 DEFINE PARAMETERS FOR THE SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the physical parameters
strength     = 0.20;           % Material parameter N0V
thouless     = 0.01;           % Thouless energy of the superconductor
upper        = 1.50;           % Upper limit on the critical temperature
lower        = 0.00;           % Lower limit on the critical temperature

% Define the simulation parameters
output       = 'output/critical_superconductor_binary/';
positions    = linspace(0, 1, 5);
energies     = [linspace(0.000,1.500,100) linspace(1.501,cosh(1/strength),100)];
iterations   = 8;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   PREPARATIONS FOR THE SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make sure that all required classes and methods are in the current path
initialize;

% Create the output directory (if it doesn't already exist)
mkdir(output);

% Log the program output to a file
diary([output, 'log.txt']);

% Create a superconductor object based on the parameters above
s = Superconductor(positions, energies, thouless, strength);

% Make sure that all debugging options are disabled
s.delay = 0;
s.debug = 0;
s.plot  = 0;

% This variable is used to keep a backup of the last non-critical object
sb = s.backup_save;

% Define variables used to store the tested temperatures and gaps
temperatures = [lower upper]; % Which temperatures we have tested
gaps         = [1 0];         % Which superconducting gaps we have found



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          PERFORM A BINARY SEARCH FOR THE CRITICAL TEMPERATURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
for n=1:iterations
    % Restore the current superconductor state from backup
    s.backup_load(sb);
    
    % Set the current temperature to the average of the two previous values
    s.temperature = (upper+lower)/2;
   
    % Status information
    fprintf(':: PROGRAM: [ %3d / %3d ] [ Temp: %.6f ] [ Time: %2d min ]\n',  n, iterations, s.temperature, floor(toc/60));
 
    % Make sure that the current gap corresponds to the state we loaded
    s.update;
    
    % Update the internal state of the superconductor until the mean gap converges
    for m=1:iterations
        if s.critical
            break;
        else
            % Update the superconductor state
            gap = s.mean_gap;
            s.update;

            % Status information
            fprintf(':: PROGRAM: [ %3d / %3d ] [ Temp: %.6f ] [ Time: %2d min ] [ Gap: %.6f ]\n            Gap changed by %.2f%%.\n',  n, iterations, s.temperature, floor(toc/60), s.mean_gap, 100*abs(1-s.mean_gap/gap));            
        end
    end
    
    % Store the current temperature and gap as results
    temperatures(end) = s.temperature;
    gaps(end)         = s.mean_gap;
    
    % If the system has gone critical, update the upper estimate for the
    % critical temperature. If not, update the lower one, and use the
    % current state as an initial guess from now on.
    if s.critical
        upper = current;
    else
        lower = current;
        sb    = s.backup_save;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           SAVE RESULTS AND CLEAN UP AFTER THE SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The final estimate of the critical temperature is the mean of the current
% upper and lower limits obtained by the above calculations
critical_temperature = (upper+lower)/2;

% Save the results to file
save([output, 'results.mat'], 'critical_temperature', 'temperatures', 'gaps');

% Disable the log from now
diary off;

% Plot the results
figure;
title('Plot of temperature vs. superconducting gap');
plot(temperatures, gaps);
xlabel('Temperature');
ylabel('Superconducting gap');