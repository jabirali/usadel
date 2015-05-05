% This script simulates the proximity effect in a normal metal connected
% to two superconductors, where there is a phase difference between them. 
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-05-06
% Updated 2015-05-06



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 DEFINE PARAMETERS FOR THE SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vectors of positions, energies, and phase differences for the simulation
positions = linspace(0,    1.0, 100);
energies  = linspace(-1.5, 1.5,  50);
phasediff = linspace(0,    pi,   50);

% Filename where results will be stored
output = 'surface_josephson.dat';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   PREPARATIONS FOR THE SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make sure that all required classes and methods are in the current path
initialize;

% Create a superconductor and normal metal
s = Superconductor([0], energies, 1, 0.2);

m = Metal(positions, energies, 1);
m.interface_left = 3;
m.interface_right = 3;
m.update_boundary_left(s);
m.update_boundary_right(s);

% This enables or disables various debugging options
m.delay = 0;
m.debug = 1;
m.plot  = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      PERFORM THE SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over phase differences, and calculate the density of states
density = zeros(length(energies),length(phasediff));
for N=1:length(phasediff)
    % Status update
    fprintf('\nPhase difference: %3.4f\n', phasediff(N));

    % Update the phase difference and state
    s.phase_set(phasediff(N));
    m.update_boundary_right(s);
    try
        m.update;
    end

    % Extract the density of states
    for M=1:length(energies)
        density(M,N) = m.states(floor(length(positions)/2), M).eval_ldos;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                POSTPROCESSING AFTER SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save the results of the simulation
save(output);

% Plot the resulting density of states as a function of polarization and phaseshift
[X,Y] = meshgrid(energies, phasediff);
surf(X,Y,density,'FaceColor','red','EdgeColor','none');
camlight('left');
lighting('phong');
