% This defines a data structure that describes the physical state of a
% metal for a given range of positions and energies. The purpose of this 
% class is mainly to be used as a base class for more interesting material
% classes, such as those that describe superconductors and ferromagnets.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-23
% Updated 2015-03-05

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
        interface_left  = inf;               % Interface parameter zeta at the left boundary
        interface_right = inf;               % Interface parameter zeta at the right boundary
        thouless        = 1;                 % Thouless energy of the system
        
        % Properties that are used during simulations
        coeff1  = {};                        % Coefficients in the differential equations for gamma
        coeff2  = {};                        % Coefficients in the differential equations for gamma~
        
        % Properties that determine the simulation behavior
        error_abs = 1e-6;                    % Maximum absolute error when simulating
        error_rel = 1e-6;                    % Maximum relative error when simulating
        grid_size = 32768;                   % Maximum grid size to use in simulations
        
        % Properties that determine the behaviour of the program
        debug         = true;                % Whether to show intermediate results or not
        plot          = true;                % Whether to plot intermediate results or not
        delay         = 0;                   % How long to wait between program iterations
    end
    
    
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods that instantiate the object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Metal(positions, energies, thouless)
            % This method constructs a Metal instance from a vector of
            % positions, a vector of energies, and the Thouless energy.
            
            % Set the internal properties to the provided values
            self.positions = positions;
            self.energies  = energies;
            self.thouless  = thouless;
            
            % Initialize the internal state to a bulk superconductor with
            % superconducting gap 1. This is useful as an initial guess
            % when simulating strong proximity effects in energy units 
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
            % the functions 'jacobian' and 'boundary' when solving equations.
            
            % Coefficients in the equations for the Riccati parameter gamma
            self.coeff1{1} = -2i/self.thouless;
            
            % Coefficients in the equations for the Riccati parameter gamma~
            self.coeff2{1} = -2i/self.thouless;
        end
        
        function update_boundary_left(self, other)
            % This function updates the boundary condition to the left
            % based on the current state of another material.
            self.boundary_left(:) = other.states(end,:);
        end

        function update_boundary_right(self, other)
            % This function updates the boundary condition to the right
            % based on the current state of another material.
            self.boundary_right(:) = other.states(1,:);
        end

        function update_state(self)
            % This function solves the Usadel equation numerically for the
            % given position and energy range, and updates the current
            % estimate for the state of the system.
            
            % Set the accuracy of the numerical solution
            options = bvpset('AbsTol',self.error_abs,'RelTol',self.error_rel,'Nmax',self.grid_size);

            % Partially evaluate the Jacobian matrix and boundary conditions
            % for the different material energies, and store the resulting 
            % anonymous functions in an array. These functions are passed on
            % to bvp6c when solving the equations.
            jc = {};
            bc = {};
            for m=1:length(self.energies)
                jc{m} = @(x,y) self.jacobian(self,x,y,self.energies(m));
                bc{m} = @(a,b) self.boundary(self,a,b,self.energies(m));
            end
                
            for m=1:length(self.energies)
                % Progress information
                self.print('[ %3.f / %3.f ]  E = %2.4f ', m, length(self.energies), self.energies(m));
                
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
                
                % Time delay between iterations (reduces load on the system)
                pause(self.delay);
            end
        end
        
        
        function update(self)
            % This function updates the internal state of the
            % metal object by calling the other update methods.
            
            % Update the state
            self.update_coeff;
            self.update_state;
                
            % Plot the current density of states (if 'plot' is set to true)
            if self.plot
                self.plot_dos;
            end
        end
    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define miscellaneous helper methods for the class
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function print(self,varargin)
            % This function is used to print messages, such as debug info
            % and progress, if the 'debug' flag is set to 'true'. The
            % message is preceded by the name of the class, which lets you
            % distinguish the output of 'Metal' instances from the output
            % of instances of daughter classes.
            
            if self.debug
                fprintf(':: %s: %s\n', class(self), sprintf(varargin{:}));
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

        function backup = backup_save(self)
            % Returns a backup of the state of the metal
            backup = self.states(:,:);
        end
        
        function backup_load(self, backup)
            % Restores the state of the metal from a backup
            % NB: Remember to run 'update' after a call to 'backup_load' to
            %     make sure that the rest of the properties are consistent!
            self.states(:,:) = backup;
        end
        
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods for printing and plotting the internal state
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function plot_dos(self)
            % Calculate the density of states for the system.
            % NB: This implementation only *adds* the contribution from
            %     every energy, and does not perform a proper integral!

            dos = zeros(length(self.energies), 1);
            for m=1:length(self.energies)
                for n=1:length(self.positions)
                    dos(m) = dos(m) + self.states(n,m).eval_ldos/length(self.positions);
                end
            end
            
            % Plot a cubic interpolation of the results
            energies = linspace(self.energies(1), self.energies(end), 100);
            plot(energies, pchip(self.energies, dos, energies));
            xlabel('Energy');
            ylabel('Density of States');
        end
        
        function plot_dist(self)
            % Calculate the singlet and triplet distributions.
            % NB: This implementation only *adds* the contribution from
            %     every energy, and does not perform a proper integral!

            singlet = zeros(length(self.positions), 1);
            triplet = zeros(length(self.positions), 1);
            for m=1:length(self.energies)
                for n=1:length(self.positions)
                    singlet(n) = singlet(n) + norm(self.states(n,m).singlet);
                    triplet(n) = triplet(n) + norm(self.states(n,m).triplet);
                end
            end
            
            % Plot cubic interpolations of the results
            positions = linspace(self.positions(1), self.positions(end), 100);
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
            
            % Extract the Riccati parameters and their derivatives
            [g,dg,gt,dgt] = State.unpack(y);
            
            % Calculate the normalization matrices
            N  = inv( eye(2) - g*gt );
            Nt = inv( eye(2) - gt*g );
            
            % Calculate the second derivatives of the Riccati parameters
            % according to the Usadel equation in the metal
            d2g  = -2 * dg*Nt*gt*dg + self.coeff1{1} * (energy+1e-3i)*g;
            d2gt = -2 * dgt*N*g*dgt + self.coeff2{1} * (energy+1e-3i)*gt;
            
            % Pack the results into a state vector
            dydx = State.pack(dg,d2g,dgt,d2gt);
        end
        
        function residue = boundary(self, y1, y2, energy)
            % This function takes a Metal object 'self', the position 'x',
            % the current state vector 'y', and an energy as inputs, and
            % calculates the Kuprianov-Lukichev boundary conditions.
            
            % State in the material to the left
            s0   = self.boundary_left(self.energy_index(energy));
            g0   = s0.g;
            dg0  = s0.dg;
            gt0  = s0.gt;
            dgt0 = s0.dgt;
            
            % State at the left end of the material
            [g1,dg1,gt1,dgt1] = State.unpack(y1);
            
            % State at the right end of the material
            [g2,dg2,gt2,dgt2] = State.unpack(y2);
            
            % State in the material to the right
            s3   = self.boundary_right(self.energy_index(energy));
            g3   = s3.g;
            dg3  = s3.dg;
            gt3  = s3.gt;
            dgt3 = s3.dgt;
            
            % Calculate the normalization matrices
            N0  = inv( eye(2) - g0*gt0 );
            Nt0 = inv( eye(2) - gt0*g0 );
            
            N3  = inv( eye(2) - g3*gt3 );
            Nt3 = inv( eye(2) - gt3*g3 );
            
            % Calculate the deviation from the Kuprianov--Lukichev b.c.
            dg1  = dg1  - (1/self.interface_left)*( eye(2) - g1*gt0 )*N0*(  g1  - g0  );
            dgt1 = dgt1 - (1/self.interface_left)*( eye(2) - gt1*g0 )*Nt0*( gt1 - gt0 );
            
            dg2  = dg2  - (1/self.interface_right)*( eye(2) - g2*gt3 )*N3*(  g2  - g3  );
            dgt2 = dgt2 - (1/self.interface_right)*( eye(2) - gt2*g3 )*Nt3*( gt2 - gt3 );
            
            % Vectorize the results of the calculations, and return it
            residue = State.pack(dg1,dgt1,dg2,dgt2);
        end
    end
end