% This script simulates the proximity effect in a metal with a spin-active interface.
% [This version of the script performs the calculation self-consistently.]
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-05-06
% Updated 2015-05-07

function simulate_bilayer_spinactive_sc(interface_polarization, interface_phase)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 DEFINE PARAMETERS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Vectors of positions and energies that will be used in the simulation
    positions     = linspace(0, 1, 150);
    energies      = [linspace(0.000,1.500,500) linspace(1.501,cosh(1/0.2),100)];
    
    % Number of iterations to perform
    iterations    = 30;
    
    % Filename where results will be stored
    output = 'simulate_bilayer_spinactive_sc.dat';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   PREPARATIONS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make sure that all required classes and methods are in the current path
    initialize;
    
    % Create a superconductor [d=3.5ξ]
    s = Superconductor(positions, energies, 1/3.5^2, 0.2);
    s.spinactive          = 1;
    s.interface_right     = 3;
    s.magnetization_right = [0,0,1];
    s.polarization_right  = interface_polarization;
    s.phaseshift_right    = interface_phase;
    
    % Create a normal metal [d=ξ]
    m = Metal(positions, energies, 1/1^2);
    m.spinactive         = 1;
    m.interface_left     = 3;
    m.magnetization_left = [0,0,1];
    m.polarization_left  = interface_polarization;
    m.phaseshift_left    = interface_phase;

    % This enables or disables various debugging options
    s.delay = 0;
    s.debug = 1;
    s.plot  = 0;
    
    m.delay = 0;
    m.debug = 1;
    m.plot  = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       PERFORM THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dosM = zeros(1,length(energies));
    dosS = zeros(1,length(energies));
    for n=1:iterations
        % Status information
        fprintf('[ %3d / %3d ] [ Time: %2d min ] [ Gap: %.6f ]\n',  n, iterations, floor(toc/60), s.gap_max);

        % Update the internal state of the normal metal
        m.update_boundary_left(s);
        m.update;
        
        % Calculate the density of states in the normal metal [left end]
        for n=1:length(energies)
            dosM(n) = m.states(1,n).eval_ldos;
        end
        disp(dosM);
    
        % Update the internal state of the superconductor
        s.update_boundary_right(m);
        s.update;

        % Calculate the density of states in the superconductor [right end]
        for n=1:length(energies)
            dosS(n) = s.states(end,n).eval_ldos;
        end
        disp(dosS);

        % Save the results of the simulation
        save(output);
    end
end