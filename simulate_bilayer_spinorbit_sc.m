% This script simulates the proximity effect in a ferromagnet with SOC.
% [This version of the script performs the calculation self-consistently.]
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-05-14
% Updated 2015-05-14

function simulate_bilayer_spinorbit_sc(exchange_strength, exchange_angle, spinorbit_strength, spinorbit_angle)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                 DEFINE PARAMETERS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Vectors of positions and energies that will be used in the simulation
    positions = linspace(0, 1, 151);
    energies  = [linspace(0.001,1.500,501) linspace(1.501,cosh(1/0.2),101)];
    
    % Filename where results will be stored
    output = 'simulate_bilayer_spinorbit_sc.dat';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   PREPARATIONS FOR THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make sure that all required classes and methods are in the current path
    initialize;

    % Create a superconductor [d=0.55ξ]
    s = Superconductor(positions, energies, 1/0.55^2, 0.2);
    s.interface_right = 3;
    %s.locked = true;
    
    % Create a ferromagnet [d=0.55ξ]
    m = Ferromagnet(positions, energies, 1/0.55^2, [exchange_strength*cos(exchange_angle), exchange_strength*sin(exchange_angle), 0], SpinVector.RashbaDresselhaus(spinorbit_strength, spinorbit_angle));
    m.interface_left = 3;
    m.update_boundary_left(s);

    % This enables or disables various debugging options
    s.delay = 0;
    s.debug = 1;
    s.plot  = 0;
    
    m.delay = 0;
    m.debug = 1;
    m.plot  = 0;
    
    % These variables keep track of the density of states 
    dosM = zeros(1,length(energies));
    dosS = zeros(1,length(energies));
    
    dosM0 = ones(1,length(energies));
    dosS0 = ones(1,length(energies));

    residue = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       PERFORM THE SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tic;
    while residue > 1e-3
        % Status information
        fprintf('[ Time: %2d min ] [ Gap: %.4f ] [ Residue: %.4f ]\n',  floor(toc/60), s.gap_max, residue);

        % Update the internal state of the normal metal
        m.update_boundary_left(s);
        m.update;
        
        % Calculate the density of states in the normal metal [left end]
        dosM0 = dosM;
        for n=1:length(energies)
            dosM(n) = m.states(1,n).eval_ldos;
        end
        fprintf('[ ZEP in ferromagnet: %.6f ]\n\n', dosM(1));
    
        % Update the internal state of the superconductor
        s.update_boundary_right(m);
        s.update;

        % Calculate the density of states in the superconductor [right end]
        dosS0 = dosS;
        for n=1:length(energies)
            dosS(n) = s.states(end,n).eval_ldos;
        end
        fprintf('[ ZEP in superconductor: %.6f ]\n\n', dosS(1));

        % Update the residue
        residue = max(max(abs(dosM-dosM0)),max(abs(dosS-dosS0)));
        
        % Save the results of the simulation
        save(output);
    end
end
