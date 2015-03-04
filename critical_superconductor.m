% This script calculates the critical temperature of a bulk superconductor,
% by performing a binary search for the temperature where the gap vanishes
% numerically. The result should be numerically one in the given unit
% system, so this script can be used to calibrate other simulation parameters.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-28
% Updated 2015-03-04


function critical_superconductor(superconductor_length, superconductor_strength)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 DEFINE PARAMETERS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Vectors of positions and energies that will be used in the simulation
    positions     = linspace(0, 1, 10);
    energies      = [linspace(0.000,1.500,500) linspace(1.501,cosh(1/superconductor_strength),100)];

    % Number of iterations of the binary search to perform
    iterations    = 8;

    % Number of iterations that the system needs to stabilize
    stabilization = 5;

    % Upper and lower limits for the binary search
    lower = 0.00;
    upper = 1.50;
    
    % Filename where results will be stored
    output = 'critical_superconductor.dat';



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   PREPARATIONS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make sure that all required classes and methods are in the current path
    initialize;

    % Create a superconductor object based on the parameters above
    s = Superconductor(positions, energies, 1/superconductor_length^2, superconductor_strength);

    % Enable or disable various debugging flags
    s.delay = 0;
    s.debug = 1;
    s.plot  = 0;

    % This variable is used to keep a backup of the last non-critical object
    sb = s.backup_save;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        PERFORM A BINARY SEARCH FOR THE CRITICAL TEMPERATURE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tic;
    for n=1:iterations
        % Set the current temperature to the average of the two previous values
        s.temperature = (upper+lower)/2;

        % Keep updating the internal state of the superconductor until we
        % either reach a phase transition, or the gap starts to increase
        loop = 0;
        gaps = [ 1 ];
        while true
            % If too many iterations have passed without convergence, then
            % proceed to accelerate the convergence by reducing the gap
            loop = loop + 1;
            if rem(loop,stabilization) == 0
                % Status information
                fprintf(':: [ %3d / %3d ] [ Temp: %.6f ] [ Time: %2d min ] Accelerating...\n',  n, iterations, s.temperature, floor(toc/60));

                % Reduce the superconducting gap everywhere in the system
                s.gap_reduce;
            end

            % Status information
            fprintf(':: [ %3d / %3d ] [ Temp: %.6f ] [ Time: %2d min ] [ Gap: %.6f ]\n',  n, iterations, s.temperature, floor(toc/60), s.gap_max);

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

            elseif (s.gap_max < 0.005)
                % The gap is so small that we must have reached critical
                % temperature. Update upper estimate, load a noncritical
                % state from backup, and terminate loop.

                upper = s.temperature;
                s.backup_load(sb);
                break;

            else
                % If the superconductor is not critical, then use the current
                % state as an initial guess in the next simulation.

                sb = s.backup_save;
            end
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           SAVE RESULTS AND CLEAN UP AFTER THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % The final estimate of the critical temperature is the mean of the current
    % upper and lower limits obtained by the above calculations
    critical_temperature = (upper+lower)/2;

    % Output the final result
    fprintf('Critical temperature: %.6f\nLower limit: %.6f\nUpper limit: %.6f\n:', critical_temperature, lower, upper);

    % Save the results to file
    save(output);
end