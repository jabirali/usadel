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
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          PERFORM A BINARY SEARCH FOR THE CRITICAL TEMPERATURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
for n=1:iterations
    % Restore the current superconductor state from backup
    s.backup_load(sb);
    
    % Set the current temperature to the average of the two previous values
    s.temperature = (upper+lower)/2;
    
    % Keep updating the internal state of the superconductor until we
    % either reach a phase transition, or the gap starts to increase
    loop = 0;
    gaps = [1];
    while true
        % If too many iterations have passed without convergence, then
        % proceed to accelerate the convergence
        loop = loop + 1;
        if rem(loop,6) == 0
            % Status information
            fprintf(':: PROGRAM: [ %3d / %3d ] [ Temp: %.6f ] [ Time: %2d min ] Accelerating...\n',  n, iterations, s.temperature, floor(toc/60));
            
            % Reduce the superconducting gap everywhere in the system
            s.gap_reduce;
        end
        
        % Status information
        fprintf(':: PROGRAM: [ %3d / %3d ] [ Temp: %.6f ] [ Time: %2d min ] [ Gap: %.6f ]\n',  n, iterations, s.temperature, floor(toc/60), s.gap_mean);
                
        % Update the superconductor state
        s.update;
        gaps(end+1) = s.gap_mean;
        
        if (gaps(end)-gaps(end-1)) > 0
            % The gap increased during the last iteration, so we must be
            % below the critical temperature. Updating the lower estimate, 
            % and using the current state as an initial guess in the next 
            % simulation. Then terminate the loop.
            
            lower = s.temperature;
            sb    = s.backup_save;
            break;

        elseif (gaps(end) < 0.01)
            % The gap is so small that we must have reached critical
            % temperature. Update upper estimate and terminate loop.

            upper = s.temperature;
            break;
            
        else
            % If the superconductor is not critical, then use the current
            % state as an initial guess in the next simulation.
            
            sb = s.backup_save;
        end
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
%figure;
%title('Plot of temperature vs. superconducting gap');
%plot(temperatures, gaps);
%xlabel('Temperature');
%ylabel('Superconducting gap');