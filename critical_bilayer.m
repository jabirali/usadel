% This script calculates the critical temperature of a superconductor/ferromagnet
% bilayer with spin-orbit coupling, by performing a kind of binary search for
% the temperature where the superconducting gap vanishes numerically. 
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-03-01
% Updated 2015-03-04



function critical_bilayer(superconductor_length, ferromagnet_length, strength, exchange, spinorbit, angle)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 DEFINE PARAMETERS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Vectors of positions and energies that will be used in the simulation
    positions     = linspace(0, 1, 100);
    energies      = [linspace(0.000,1.500,500) linspace(1.501,cosh(1/strength),100)];

    % Number of iterations of the binary search to perform
    iterations    = 8;
    
    % Number of iterations that the system needs to stabilize
    stabilization = 8;

    % Upper and lower limits for the binary search
    lower = 0;
    upper = 1;

    % Filename where results will be stored
    output = 'critical_bilayer.dat';



    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   PREPARATIONS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make sure that all required classes and methods are in the current path
    initialize;
    
    % Instantiate and initialize superconductor/ferromagnet objects
    s = Superconductor(positions, energies, 1/superconductor_length^2, strength);
    f = Ferromagnet(positions, energies, 1/ferromagnet_length^2, exchange, SpinVector.RashbaDresselhaus(spinorbit, angle));
    
    s.temperature     = 0;
    s.interface_right = 3;
    f.interface_left  = 3;

    % This enables or disables various debugging options
    s.delay = 0;
    s.debug = 1;
    s.plot  = 0;
    f.delay = 0;
    f.debug = 1;
    f.plot  = 0;
    
    % Print out parameters for verification purposes
    fprintf('SUPERCONDUCTOR:\n :: Length:          %1.6f\n :: Thouless energy: %1.6f\n :: Interface:       %1.6f\n :: Strength:        %1.6f\n\n', superconductor_length, s.thouless, s.interface_right, s.strength);
    fprintf('FERROMAGNET:\n :: Length:          %1.6f\n :: Thouless energy: %1.6f\n :: Interface:       %1.6f\n :: Exchange field h:\n', ferromagnet_length, f.thouless, f.interface_left);
    disp(f.exchange);
    fprintf(' :: Spin-orbit field Ax:\n')
    disp(f.spinorbit.x);
    fprintf(' :: Spin-orbit field Ay:\n')
    disp(f.spinorbit.y);
    fprintf(' :: Spin-orbit field Az:\n')
    disp(f.spinorbit.z);

    
    % Initialize the bilayer by performing 'stabilization' iterations at
    % zero temperature, to make sure that we do get the proximity effect
    tic;
    for n=1:stabilization
        % Status information
        fprintf('[ %3d / %3d ] [ Temp: %.6f ] [ Time: %2d min ] Initializing...\n',  n, stabilization, s.temperature, floor(toc/60));
        
        % Update the boundary condition and state of the ferromagnet
        f.update_boundary_left(s);
        f.update;
        
        % Update the boundary condition and state of the superconductor
        s.update_boundary_right(f);
        s.update;
    end

    % These variables are used to keep a backup of the last non-critical object
    sb = s.backup_save;
    fb = f.backup_save;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %          PERFORM A BINARY SEARCH FOR THE CRITICAL TEMPERATURE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
                fprintf('[ %3d / %3d ] [ Temp: %.6f ] [ Time: %2d min ] Accelerating...\n',  n, iterations, s.temperature, floor(toc/60));

                % Reduce the superconducting gap everywhere in the system
                s.gap_reduce;
            end

            % Status information
            fprintf('[ %3d / %3d ] [ Temp: %.6f ] [ Time: %2d min ] [ Gap: %.6f ]\n',  n, iterations, s.temperature, floor(toc/60), s.gap_max);

            % Update the ferromagnet boundary conditions and state
            f.update_boundary_left(s);
            f.update;
            
            % Update the superconductor boundary conditions and state
            s.update_boundary_right(f);
            s.update;
            
            % Store the current superconductor mean gap in 'gaps' 
            gaps(end+1) = s.gap_mean;

            % This is the logic that controls the while loop
            if (gaps(end)-gaps(end-1)) > 0
                % The gap increased during the last iteration, so we must be
                % below the critical temperature. Updating the lower estimate, 
                % and using the current state as an initial guess in the next 
                % simulation. Then terminate the loop.

                lower = s.temperature;
                sb    = s.backup_save;
                fb    = f.backup_save;
                break;

            elseif (s.gap_max < 0.005)
                % The gap is so small that we must have reached critical
                % temperature. Update upper estimate, load a noncritical
                % state from backup, and terminate the loop.

                upper = s.temperature;
                s.backup_load(sb);
                f.backup_load(fb);
                break;

            else
                % If the superconductor is not critical, then use the current
                % state as an initial guess in future simulations.

                sb = s.backup_save;
                fb = f.backup_save;
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