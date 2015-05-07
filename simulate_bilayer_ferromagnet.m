% This script simulates the proximity effect in a ferromagnet with SOC.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-06
% Updated 2015-05-07

function simulate_bilayer_ferromagnet(exchange_strength, exchange_angle, spinorbit_strength, spinorbit_angle)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 DEFINE PARAMETERS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Vectors of positions and energies that will be used in the simulation
    positions     = linspace(0.0, 1.0, 100);
    energies      = linspace(0.0, 1.5,  25);
    
    % Filename where results will be stored
    output = 'simulate_bilayer_ferromagnet.dat';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   PREPARATIONS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make sure that all required classes and methods are in the current path
    initialize;
    
    % Create a superconductor
    s = Superconductor([0], energies, 1, 0.2);
    
    % Create a ferromagnet connected to the superconductor above
    m = Ferromagnet(positions, energies, 1, [exchange_strength*cos(exchange_angle), exchange_strength*sin(exchange_angle), 0], SpinVector.RashbaDresselhaus(spinorbit_strength, spinorbit_angle));
    m.interface_left = 3;
    m.update_boundary_left(s);
    
    % This enables or disables various debugging options
    m.delay = 0;
    m.debug = 1;
    m.plot  = 0;

    % Update the internal state of the metal
    m.update;
    
    % Plot the results
    figure;
    m.plot_dos_surf;
    
    % Save the results of the simulation
    save(output);
end
