% This script simulates the proximity effect in a metal with a spin-active interface.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-05-06
% Updated 2015-05-07

function simulate_bilayer_spinactive(interface_polarization, interface_phase)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 DEFINE PARAMETERS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Vectors of positions and energies that will be used in the simulation
    positions     = linspace(0.0, 1.0, 100);
    energies      = linspace(0.0, 1.5,  26);
    
    % Filename where results will be stored
    output = 'simulate_bilayer_spinactive.dat';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   PREPARATIONS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make sure that all required classes and methods are in the current path
    initialize;
    
    % Create a superconductor
    s = Superconductor([0], energies, 1, 0.2);

    % Create a normal metal connected to the superconductors above
    m = Metal(positions, energies, 1);
    m.spinactive         = 1;
    m.interface_left     = 3;
    m.magnetization_left = [0,0,1];
    m.polarization_left  = interface_polarization;
    m.phaseshift_left    = interface_phase;
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