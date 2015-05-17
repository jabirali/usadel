% This script simulates the proximity effect in a normal metal connected
% to two superconductors with a constant phase difference between them.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-05-06
% Updated 2015-05-07

function simulate_josephson(phase_difference)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 DEFINE PARAMETERS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Vectors of positions and energies that will be used in the simulation
    positions     = linspace(0.0, 1.0, 50);
    energies      = linspace(0.0, 2.0, 50);
    
    % Filename where results will be stored
    output = 'simulate_josephson.dat';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   PREPARATIONS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make sure that all required classes and methods are in the current path
    initialize;
    
    % Create two superconductors with 'phase_difference' between them
    s1 = Superconductor([0], energies, 1, 0.2);
    s1.complex = true;
    s1.phase_set(-phase_difference/2);
    
    s2 = Superconductor([0], energies, 1, 0.2);
    s2.complex = true;
    s2.phase_set(+phase_difference/2);

    % Create a normal metal connected to the superconductors above
    m = Metal(positions, energies, 1);
    m.interface_left  = 3;
    m.interface_right = 3;
    m.update_boundary_left(s1);
    m.update_boundary_right(s2);
    
    % This enables or disables various debugging options
    m.delay = 0;
    m.debug = 1;
    m.plot  = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       PERFORM THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Update the internal state of the metal
    m.update;
    
    % Plot the results
    figure;
    m.plot_dos_surf;
    figure;
    m.plot_dos_center;
    
    % Save the results of the simulation
    save(output);
end