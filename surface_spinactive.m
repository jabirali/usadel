% This script simulates the proximity effect in a normal metal connected
% to two superconductors, where one interface is assumed to be spin-active. 
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-05-05
% Updated 2015-05-06



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 DEFINE PARAMETERS FOR THE SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vectors of positions that will be used in the simulation
positions     = linspace(0, 1.0, 200);

% Filename where results will be stored
output = 'surface_spinactive.dat';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   PREPARATIONS FOR THE SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make sure that all required classes and methods are in the current path
initialize;

% Create a metal with one spin-active interface and one regular interface
m = Metal(positions, [0], 1/3^2);
m.spinactive = true;

m.interface_left      = 3;
m.magnetization_left  = [0,0,1];
m.boundary_left(1)    = Superconductor.Bulk(0,1,0);

m.interface_right     = 3;
m.magnetization_right = [0,0,0];
m.boundary_right(1)   = Superconductor.Bulk(0,1,0);

% This enables or disables various debugging options
m.delay = 0;
m.debug = 0;
m.plot  = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      PERFORM THE SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over polarizations and phaseshifts, and calculate the zero-energy
% density of states at the superconductor/normal-metal interface
density      = zeros(20,20);
polarization = linspace(-1.0,1.0,20);
phaseshift   = linspace(-1.0,1.0,20);

for N=1:20
    m.polarization_left = polarization(N);
    for M=1:20
        m.phaseshift_left = phaseshift(M);
        fprintf('[ Polarization: %3.4f ] [ Phaseshift: %3.4f ]\n', polarization(N), phaseshift(M));
        try
            m.update;
        end
        
        density(N,M) = m.states(1,1).eval_ldos;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                POSTPROCESSING AFTER SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the resulting density of states as a function of polarization and phaseshift
[X,Y] = meshgrid(polarization, phaseshift);
surf(X,Y,density,'FaceColor','red','EdgeColor','none');
camlight('left');
lighting('phong');

% Save the results of the simulation
save(output);