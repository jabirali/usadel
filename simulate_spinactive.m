% This script simulates the proximity effect in a bilayer that consists of
% a superconductor and normal metal connected by a spin-active interface.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-05-05
% Updated 2015-05-05

function simulate_spinactive(superconductor_length, metal_length, polarization, phaseshift)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 DEFINE PARAMETERS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Vectors of positions and energies that will be used in the simulation
    positions     = linspace(0, 1.0, 100);
    energies      = linspace(0, 1.5, 100);

    % Filename where results will be stored
    output = 'simulate_spinactive.dat';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   PREPARATIONS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make sure that all required classes and methods are in the current path
    initialize;

    % Create a superconductor and normal metal based on the parameters above
    s = Superconductor(positions, energies, 1/superconductor_length^2, 0.2);
    m = Metal(positions, energies, 1/metal_length^2);
    
    % Use a spin-active interface as the inner boundary conditions
    s.spinactive          = true;
    s.interface_right     = 3;
    s.polarization_right  = polarization;
    s.phaseshift_right    = phaseshift;
    s.magnetization_right = [0,0,1];
    s.update_boundary_right(m);

    m.spinactive         = true;
    m.interface_left     = 3;
    m.polarization_left  = polarization;
    m.phaseshift_left    = phaseshift;
    m.magnetization_left = [0,0,1];
    m.update_boundary_left(s);
    
    % Use bulk superconductors as the outer boundary conditions
    s.interface_left = 3;
    s.update_boundary_left(s);
    
    m.interface_right = 3;
    m.update_boundary_right(s);
    
    % This enables or disables various debugging options
    m.delay = 0;
    m.debug = 1;
    m.plot  = 1;
    
    s.delay = 0;
    s.debug = 1;
    s.plot  = 1;

    % Update the internal state of the materials
    m.update;

    % Postpone self-consistent calculation until parameters determined
%     for n=1:8
%         s.update_boundary_right(m);
%         s.update;
%     
%         m.update_boundary_left(s);
%         m.update;
%     end
    
    % Plot the results
    figure;
    s.plot_dos;
    figure;
    s.plot_dist;

    figure;
    m.plot_dos;
    figure;
    m.plot_dist;
    
    % Save the results of the simulation
    save(output);
end