% This script simulates the proximity effect in a ferromagnet, where the
% left end of the ferromagnet is connected to a bulk superconductor.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-03-04
% Updated 2015-03-04

function simulate_ferromagnet(ferromagnet_length, exchange, spinorbit, angle)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 DEFINE PARAMETERS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Vectors of positions and energies that will be used in the simulation
    positions     = linspace(0, 1, 32);
    energies      = linspace(0, 2, 100);

    % Filename where results will be stored
    output = 'simulate_ferromagnet.dat';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   PREPARATIONS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make sure that all required classes and methods are in the current path
    initialize;

    % Create a ferromagnet based on the parameters above
    f = Ferromagnet(positions, energies, 1/ferromagnet_length^2, exchange, SpinVector.RashbaDresselhaus(spinorbit, angle));

    % This enables or disables various debugging options
    f.delay = 0;
    f.debug = 1;
    f.plot  = 0;

    % Set the boundary conditions for the ferromagnet (i.e. connect it to a bulk superconductor)
    %f.transparent = true;
    f.interface_left  = 1;
    f.interface_right = inf;
    for m=1:length(f.energies)
        f.boundary_left(m) = Superconductor.Bulk(f.energies(m), 1.0, 0);
    end

    % Update the internal state of the ferromagnet
    f.update;
    
    % Plot the results
    figure;
    f.plot_dos;
    figure;
    f.plot_dist;
    
    % Save the results of the simulation
    save(output);
end