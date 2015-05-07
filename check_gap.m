% This script calculates the value of the superconducting gap for an
% SF bilayer with spin-orbit coupling at some given temperature.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-05-05
% Updated 2015-05-06



function check_gap(superconductor_length, ferromagnet_length, strength, exchange, spinorbit, angle, temperature)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 DEFINE PARAMETERS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Vectors of positions and energies that will be used in the simulation
    positions  = linspace(0, 1, 150);
    energies   = [linspace(0.000,1.500,500) linspace(1.501,cosh(1/strength),100)];

    % Number of iterations that the system needs to stabilize
    iterations = 50;
    


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   PREPARATIONS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make sure that all required classes and methods are in the current path
    initialize;
    
    % Instantiate and initialize superconductor/ferromagnet objects
    s = Superconductor(positions, energies, 1/superconductor_length^2, strength);
    f = Ferromagnet(positions, energies, 1/ferromagnet_length^2, exchange, SpinVector.RashbaDresselhaus(spinorbit, angle));
    
    s.temperature     = temperature;
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

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                      PERFORM THE CALCULATIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Perform a set number of iterations at the given temperature to
    % simulate the proximity effect in the bilayer self-consistently
    tic;
    for n=1:iterations
        % Status information
        fprintf('[ %3d / %3d ] [ Temp: %.6f ] [ Time: %2d min ] [ Gap: %.6f ]\n',  n, iterations, s.temperature, floor(toc/60), s.gap_max);
        
        % Update the boundary condition and state of the ferromagnet
        f.update_boundary_left(s);
        f.update;
        
        % Update the boundary condition and state of the superconductor
        s.update_boundary_right(f);
        s.update;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                         POSTPROCESSING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Output the final result
    fprintf('Temperature: %.6f\nGap: %.6f\n', s.temperature, s.gap_max);
end
