% This script calculates the critical temperature of a superconductor/ferromagnet
% bilayer with spin-orbit coupling, by performing a binary search for the 
% temperature where the superconducting gap vanishes numerically. 
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-03-01
% Updated 2015-03-02



function critical_bilayer_binary(thouless, strength, exchange, spinorbit)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 DEFINE PARAMETERS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define the simulation parameters
    output        = 'output/critical_bilayer/';
    positions     = linspace(0, 1, 5);
    energies      = [linspace(0.000,1.500,500) linspace(1.501,cosh(1/strength),100)];
    iterations    = 8;  % Number of iterations of the binary search to perform
    stabilization = 8;  % Number of iterations the system needs to stabilize


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   PREPARATIONS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make sure that all required classes and methods are in the current path
    initialize;

    % Create the output directory (if it doesn't already exist)
    mkdir(output);

    % Log the program output to a file
    diary([output, 'log.txt']);

    % Instantiate and initialize superconductor/ferromagnet objects
    s = Superconductor(positions, energies, thouless, strength);
    f = Ferromagnet(positions, energies, thouless, exchange, spinorbit);
    
    s.temperature     = 0;
    s.interface_right = 3;
    f.interface_left  = 3;

    % Make sure that all debugging options are disabled
    s.delay = 0;
    s.debug = 0;
    s.plot  = 0;
    f.delay = 0;
    f.debug = 0;
    f.plot  = 0;
    
    % Initialize the bilayer by performing 'stabilization' iterations at
    % zero temperature, to make sure that we get a proximity effect
    for n=1:stabilization
        % Status information
        fprintf('[ %3d / %3d ] [ Temp: %2d min ] [ Time: %2d min ] Initializing state...\n',  n, stabilization, s.temperature, floor(toc/60));
        
        % Update the boundary condition and state of the ferromagnet
        f.update_boundary_left(s);
        f.update;
        
        % Update the boundary condition and state of the superconductor
        s.update_boundary_right(f);
        s.update;
    end


    % This variable is used to keep a backup of the last non-critical object
    sb = s.backup_save;
    fb = f.backup_save;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %          PERFORM A BINARY SEARCH FOR THE CRITICAL TEMPERATURE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tic;
    for n=1:iterations
        % Restore the current state from backup
        s.backup_load(sb);
        f.backup_load(fb);

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
            if rem(loop,stabilization) == 0
                % Status information
                fprintf('[ %3d / %3d ] [ Temp: %.6f ] [ Time: %2d min ] Accelerating...\n',  n, iterations, s.temperature, floor(toc/60));

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
                fb    = f.backup_save;
                break;

            elseif (gaps(end) < 0.005)
                % The gap is so small that we must have reached critical
                % temperature. Update upper estimate and terminate loop.

                upper = s.temperature;
                break;

            else
                % If the superconductor is not critical, then use the current
                % state as an initial guess in the next simulation.

                sb = s.backup_save;
                fb = f.backup_save;
            end
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           SAVE RESULTS AND CLEAN UP AFTER THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % The final estimate of the critical temperature is the mean of the current
    % upper and lower limits obtained by the above calculations
    critical_temperature = (upper+lower)/2;

    % Output the final result
    fprintf('Critical temperature: %.6f\nLower limit: %.6f\nUpper limit: %.6f\n:', critical_temperature, lower, upper);

    % Save the results to file
    save([output, 'result.mat'], 'critical_temperature', 'upper', 'lower');

    % Disable the log from now
    diary off;
end