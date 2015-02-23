% This defines a data structure that describes the physical state of a
% metal for a given range of positions and energies. The most demanding
% calculations are parallellized using SPMD. The purpose of this class is
% mainly to be used as a base class for more interesting material classes,
% such as superconductors and ferromagnets.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-23
% Updated 2015-02-23

classdef Metal < handle
    properties (GetAccess=public, SetAccess=public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define the internal variables of the Metal class
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Properties used to store the internal state of the system
        positions       = [];                % Vector of positions in the material
        energies        = [];                % Vector of energies for the material
        states          = State.empty(0,0);  % Riccati parameters and derivatives for each position and energy
        boundary_left   = State.empty(0);    % Riccati parameters and derivatives at the left boundary
        boundary_right  = State.empty(0);    % Riccati parameters and derivatives at the right boundary
        
        % Properties that determine the physical characteristics of the system
        length          = 1;                 % Length of the system
        thouless        = 1;                 % Thouless energy of the system
        interface_left  = inf;               % Interface parameter zeta at the left boundary
        interface_right = inf;               % Interface parameter zeta at the right boundary
        
        % Properties that are used during simulations
        coeff1  = {};                        % Coefficients in the differential equations for gamma
        coeff2  = {};                        % Coefficients in the differential equations for gamma~
        
        % Properties that determine the simulation behavior
        sim_error_abs = 1e-3;                % Maximum absolute error when simulating
        sim_error_rel = 1e-3;                % Maximum relative error when simulating
        sim_grid_size = 512;                 % Maximum grid size to use in simulations
        
        % Properties that determine the behaviour of the program
        debug         = true;                % Whether to show intermediate results or not
        plot          = true;                % Whether to plot intermediate results or not
        delay         = 0;                   % How long to wait between program iterations
    end
    
    
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods that instantiate the object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Metal(positions, energies, material_length, material_thouless)
            % This method constructs a Metal instance from a vector of
            % positions, a vector of energies, the length of the material,
            % and the Thouless energy of the material.
            
            % Set the internal properties to the provided values
            self.positions = positions;
            self.energies  = energies;
            self.length    = material_length;
            self.thouless  = material_thouless;
            
            % Initialize the internal state to a bulk superconductor with
            % superconducting gap 1. This is useful as an initial guess
            % when simulating strong proximity effects with energy units 
            % where the zero-temperature gap is normalized to unity.
            self.states(length(positions), length(energies)) = State;
            for i=1:length(positions)
                for j=1:length(energies)
                    self.states(i,j) = Superconductor.Bulk(energies(j), 1);
                end
            end
            
            % Set the boundary conditions to vacuum states
            self.boundary_left(length(energies))  = State;
            self.boundary_right(length(energies)) = State;
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods that update the internal state of the object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function update_coeff(self)
            % This function updates the vector of coefficients passed to
            % the differential equation solver, which in turn passes it to
            % the functions 'usadel' and 'boundary' when solving equations.
            
            % Coefficients in the equations for the Riccati parameter gamma
            self.coeff1{1} = -2;
            self.coeff1{2} = -2i/self.thouless;
            
            % Coefficients in the equations for the Riccati parameter gamma~
            self.coeff2{1} = self.coeff1{1};
            self.coeff2{2} = self.coeff1{2};
        end

        function update_state(self)
            % This function solves the Usadel equation numerically for the
            % given position and energy range, and updates the current
            % stored state of the system.
            
            % Set the accuracy of the numerical solution
            options = bvpset('AbsTol',self.sim_error_abs,'RelTol',self.sim_error_rel,'Nmax',self.sim_grid_size);

            % Partially evaluate the Jacobian matrix and boundary conditions
            % for the different material energies, and store the resulting 
            % anonymous functions in a vector. These functions are passed on
            % to bvp6c when solving the equations.
            jc = {};
            bc = {};
            for m=1:length(self.energies)
                jc{m} = @(x,y) self.jacobian(self,x,y,self.energies(m));
                bc{m} = @(a,b) self.boundary(self,a,b,self.energies(m));
            end
                
            % Parallelize the loop over the energies of the system
            %spmd
                for m=drange(1:length(self.energies))
                    % Progress information
                    self.print('[ %2.f / %2.f ]   iteration starting...', m, length(self.energies));
                    
                    % Vectorize the current state of the system for the given
                    % energy, and use it as an initial guess for the solution
                    current = zeros(16,length(self.positions));
                    for n=1:length(self.positions)
                        current(:,n) = self.states(n,m).vectorize;
                    end
                    initial = bvpinit(self.positions', current);
                    
                    % Solve the differential equation, and evaluate the
                    % solution on the position vector of the metal
                    solution = deval(bvp6c(jc{m},bc{m},initial,options), self.positions);

                    % Update the current state of the system based on the solution
                    for n=1:length(self.positions)
                        self.states(n,m) = State(solution(:,n));
                    end
                                    
                    % Progress information
                    self.print('[ %2.f / %2.f ]   iteration complete!', m, length(self.energies));
                    
                    % Small time delay to prevent the interpreter from getting sluggish or killed by the system
                    pause(self.delay);
                    %end
                end
            end
        
        
        function update(self)
            % This function updates the internal state of the
            % metal object by calling the other update methods.
            
            % Update the state
            self.update_coeff;
            self.update_state;
                
            % Plot the current DOS
            if self.plot
                self.plot_dos;
            end
        end
    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define miscellaneous helper methods for the class
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function print(self,varargin)
            % This function is used to print messages, such as debug info
            % and progress, if the 'debug' flag is set to 'true'.
            
            if self.debug
                fprintf(':: %s\n', sprintf(varargin{:}));
            end
        end

        function index = position_index(self, position)
            % Returns the vector index corresponding to a given energy value
            index = find(abs(self.positions-position) < 1e-8, 1, 'first');
        end
        
        function index = energy_index(self, energy)
            % Returns the vector index corresponding to a given energy value
            index = find(abs(self.energies-energy) < 1e-8, 1, 'first');
        end

        
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods for printing and plotting the internal state
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function plot_dos(self)
            % Calculate the density of states
            dos = zeros(length(self.energies), 1);
            for m=1:length(self.energies)
                for n=1:length(self.positions)
                    dos(m) = dos(m) + self.states(n,m).eval_ldos;
                end
            end
                        
            % Plot a cubic interpolation of the results
            energies = linspace(0,self.energies(end), 100);
            plot(energies, pchip(self.energies, dos, energies));
            xlabel('Energy');
            ylabel('Density of States');
        end
        
        function plot_dist(self)
            % Calculate the singlet and triplet distributions
            singlet = zeros(length(self.positions), 1);
            triplet = zeros(length(self.positions), 1);
            for m=1:length(self.energies)
                for n=1:length(self.positions)
                    singlet(n) = singlet(n) + norm(self.states(n,m).singlet);
                    triplet(n) = triplet(n) + norm(self.states(n,m).triplet);
                end
            end
                        
            % Plot cubic interpolations of the results
            positions = linspace(0, self.positions(end), 100);
            plot(positions, pchip(self.positions, singlet, positions), ...
                 positions, pchip(self.positions, triplet, positions));
            xlabel('Relative position');
            ylabel('Distribution');
            legend('Singlet', 'Triplet');
        end
    end
    
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define static methods (available without object instantiation)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function dydx = jacobian(self, x, y, energy)
            % This function takes a Metal object 'self', the position 'x', 
            % the current state vector 'y', and an energy as inputs, and 
            % calculates the Jacobian of the system. This is performed 
            % using the Riccati parametrized Usadel equations.
            
            % Instantiate a 'State' object based on the state vector
            state = State(y);
                                    
            % Extract the Riccati parameters and their derivatives
            g   = state.g;
            dg  = state.dg;
            gt  = state.gt;
            dgt = state.dgt;
            
            % Calculate the normalization matrices
            N  = inv( eye(2) - g*gt );
            Nt = inv( eye(2) - gt*g );
            
            % Calculate the second derivatives of the Riccati parameters
            % according to the Usadel equation in the metal
            d2g  = self.coeff1{1} * dg*Nt*gt*dg + self.coeff1{2} * energy*g;
            d2gt = self.coeff2{1} * dgt*N*g*dgt + self.coeff2{2} * energy*gt;
                        
            % Fill the results of the calculations back into a 'State' object
            state.g   = dg;
            state.dg  = d2g;
            state.gt  = dgt;
            state.dgt = d2gt;
            
            % Pack the results into a state vector
            dydx = state.vectorize;
        end
        
        function residue = boundary(self, y1, y2, energy)
            % This function takes a Metal object 'self', the position 'x', 
            % the current state vector 'y', and an energy as inputs, and 
            % calculates the Kuprianov-Lukichev boundary conditions.
            
            % State in the material to the left
            s0   = self.boundary_left(self.energy_index(energy));
            
            % State at the left end of the material
            s1   = State(y1);
            
            % State at the right end of the material
            s2   = State(y2);
            
            % State in the material to the right
            s3   = self.boundary_right(self.energy_index(energy));
             
            % Calculate the normalization matrices
            N0  = inv( eye(2) - s0.g*s0.gt );
            Nt0 = inv( eye(2) - s0.gt*s0.g );

            N3  = inv( eye(2) - s3.g*s3.gt );
            Nt3 = inv( eye(2) - s3.gt*s3.g );
            
            % Calculate the deviation from the Kuprianov--Lukichev boundary
            % conditions, and store the results back into State instances
            s1.dg  = s1.dg  - ( eye(2) - s1.g*s0.gt )*N0*(  s1.g  - s0.g  )/self.interface_left;
            s1.dgt = s1.dgt - ( eye(2) - s1.gt*s0.g )*Nt0*( s1.gt - s0.gt )/self.interface_left;
            
            s2.dg  = s2.dg  - ( eye(2) - s2.g*s3.gt )*N3*(  s2.g  - s3.g  )/self.interface_right;
            s2.dgt = s2.dgt - ( eye(2) - s2.gt*s3.g )*Nt3*( s2.gt - s3.gt )/self.interface_right;

            % Vectorize the results of the calculations, and return it            
            residue = [s1.vectorize_dg s1.vectorize_dgt s2.vectorize_dg s2.vectorize_dgt]';
        end
    end
end