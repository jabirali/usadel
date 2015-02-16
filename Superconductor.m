% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-15
% Updated 2015-02-16
%
% This defines a data structure that describes the physical state of
% superconducting material for a given range of positions and energies.


classdef Superconductor < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the internal variables for the data structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (GetAccess=public, SetAccess=public)
        positions   = [];                    % Positions in the superconductor
        energies    = [];                    % Energies of the superconductor
        gap         = [];                    % Superconducting gap at each position
        states      = State.empty(0,0);      % Green's functions for each position and energy
        
        temperature     = 1e-16;             % Temperature of the system 
        scaling         = 1;                 % Material constant N₀λ
        diffusion       = 1;                 % Diffusion constant
        interface_left  = inf;               % Interface parameter zeta (left)
        interface_right = inf;               % Interface parameter zeta (right)
        
        boundary_left   = State.empty(0);    % Boundary condition (left)   (default: vacuum)
        boundary_right  = State.empty(0);    % Boundary condition (right)  (default: vacuum)
        
        sim_error_abs = 1e-2;                % Maximum absolute error when simulating
        sim_error_rel = 1e-2;                % Maximum relative error when simulating
        sim_grid_size = 256;                 % Maximum grid size to use in simulations
    end
    

    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define constructor and accessor methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

        function self = Superconductor(positions, energies)
            % Define a constructor which initializes the Superconductor
            % from a vector of positions and a vector of energies
            self.positions = positions;
            self.energies  = energies;
            self.gap       = ones(size(positions));
            
            % Initialize the internal state to a BCS bulk superconductor
            self.states(length(positions), length(energies)) = 0;
            for i=1:length(positions)
                for j=1:length(energies)
                    self.states(i,j) = Superconductor.Bulk(energies(j), 1);
                end
            end
            
            % Set the boundary conditions to vacuum states by default
            self.boundary_left(length(energies))  = 0;
            self.boundary_right(length(energies)) = 0;    
        end
        
        function index = position_index(self, position)
            % Returns the vector index corresponding to a given energy value
            index = find(abs(self.positions-position) < 1e-5, 1, 'first');
            if ~isscalar(index)
                error('Superconductor.position_index: Provided value is not in the position vector!');
            end
        end

        function index = energy_index(self, energy)
            % Returns the vector index corresponding to a given energy value
            index = find(abs(self.energies-energy) < 1e-5, 1, 'first');
            if ~isscalar(index)
                error('Superconductor.energy_index: Provided value is not in the energy vector!');
            end
        end
                
        function result = critical(self)
            % This function returns whether or not the system is above
            % critical temperature, i.e. if the superconducting gap is zero
            % everywhere in the superconductor.
            result = ( max(abs(self.gap)) < 1e-4 );
        end
        
        function result = gap_lookup(self, x)
            % This function performs a nearest-neighbor interpolation of
            % the superconducting gap as a function of position, and
            % returns the value in a given point.
            
            result = interp1(self.positions, self.gap, x, 'nearest', 'extrap');
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods that update the internal state of the object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        
        function gap_update(self)
            % This function updates the vector containing the current
            % estimate of the superconducting gap at equilibrium.
            
            for n=1:length(self.positions)
                self.gap(n) = self.gap_calculate(self, self.positions(n));
            end
        end
        
        function state_update(self)
            % This function solves the Usadel equation numerically for the
            % given position and energy range, using the current stored 
            % estimate for the superconducting gap.

            % Set the accuracy of the numerical solution
            options = bvpset('AbsTol',self.sim_error_abs,'RelTol',self.sim_error_rel,'Nmax',self.sim_grid_size);
            
            % Information about parallel execution
            task = getCurrentTask();
            if isempty(task)
                taskID = 0;
            else
                taskID = task.ID;
            end

            % Start a timer to predict ETA
            tic;

            for m=1:length(self.energies)
                % Vectorize the current state of the system for the given
                % energy, and use it as an initial guess for the solution
                current = zeros(16,length(self.positions));
                for n=1:length(self.positions)
                    current(:,n) = self.states(n,m).vectorize;
                end
                initial = bvpinit(self.positions', current);
                
                % Partially evaluate the Jacobian and boundary conditions
                % for the current superconductor energy
                jc = @(x,y) Superconductor.jacobian(self,x,y,self.energies(m));
                bc = @(a,b) Superconductor.boundary(self,a,b,self.energies(m));
                
                % Solve the differential equation, and evaluate the
                % solution on the position vector of the superconductor 
                solution = deval(bvp6c(jc,bc,initial,options), self.positions);
                
                % Update the current state of the system based on the solution
                for n=1:length(self.positions)
                    self.states(n,m) = State(solution(:,n));
                end
                
                % Progress information
                disp(sprintf('-- Worker %2.0f: ETA: %2.f sec [ %2.f / %2.f ]', taskID, toc*(1-m/length(self.energies)), m, length(self.energies)));
                pause(0.1);
            end
        end
        
        function update(self)
            % This function updates the internal state of the
            % superconductor object by calling first 'state_update' and 
            % then 'gap_update'. Always run this after changing the
            % boundary conditions or temperature of the system.
            
            self.state_update;
            self.gap_update;
        end
    end 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define static methods (available without object instantiation)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    methods (Static)        
        function result = Bulk(energy, gap)
            % This function takes as its input an energy and a superconducting gap,
            % and returns a State object with Green's functions that correspond to
            % a BCS superconductor bulk state.
            
            theta = atanh(gap/(energy+0.001i));
            c     = cosh(theta);
            s     = sinh(theta);
            
            result = State([0,  s/(1+c); -s/(1+c), 0], 0,            ...
                           [0, -s/(1+c);  s/(1+c), 0], 0);
        end
    
        function dydx = jacobian(self, x, y, energy)
            % This function takes a Superconductor object 'self', the
            % position 'x', the current state vector 'y', and an energy as
            % inputs, and calculates the Jacobian of the system. This is
            % performed using the Riccati parametrized Usadel equations.
            
            % Instantiate a 'State' object based on the state vector
            state = State(y);
            
            % Extract diffusion constant and superconducting gap
            diff = self.diffusion;
            gap  = self.gap_lookup(x);
            
            % Extract the Riccati parameters and their derivatives
            g   = state.g;
            dg  = state.dg;
            gt  = state.gt;
            dgt = state.dgt;
            
            % Calculate the normalization matrices
            N  = inv( eye(2) - g*gt );
            Nt = inv( eye(2) - gt*g );
            
            % Calculate the second derivatives of the Riccati parameters
            % according to the Usadel equation in the superconductor
            d2g  =  - 2*dg*Nt*gt*dg ...
                    - 2i*((energy+0.001i)/diff)*g   ...
                    - (gap/diff)*(SpinVector.Pauli.y - g * SpinVector.Pauli.y * g);
            
            d2gt =  - 2*dgt*N*g*dgt  ...
                    - 2i*((energy+0.001i)/diff)*gt   ...
                    + (gap/diff)*(SpinVector.Pauli.y - gt * SpinVector.Pauli.y * gt);
            
            % Fill the results of the calculations back into a 'State' object
            state.g   = dg   + 0.0001i;
            state.dg  = d2g  + 0.0001i;
            state.gt  = dgt  - 0.0001i;
            state.dgt = d2gt - 0.0001i;
            
            % Pack the results into a state vector
            dydx = state.vectorize;
        end
        
        function residue = boundary(self, y1, y2, energy)
            % This function takes a Superconductor object 'self', the
            % position 'x', the current state vector 'y', and an energy as
            % inputs, and calculates the Kuprianov-Lukichev boundary
            % conditions for the system. 
            
            % State in the material to the left of the superconductor
            s0   = State(self.boundary_left(self.energy_index(energy)));
            
            % State at the left end of the superconductor
            s1   = State(y1);
            
            % State at the right end of the superconductor
            s2   = State(y2);
            
            % State in the material to the right of the superconductor
            s3   = State(self.boundary_right(self.energy_index(energy)));
             
            % Calculate the normalization matrices
            N0  = inv( eye(2) - s0.g*s0.gt );
            Nt0 = inv( eye(2) - s0.gt*s0.g );

            N3  = inv( eye(2) - s3.g*s3.gt );
            Nt3 = inv( eye(2) - s3.gt*s3.g );
            
            % Calculate interface parameters in the Kuprianov-Lukichev B.C.
            param_left  = abs(self.positions(end) - self.positions(1)) * self.interface_left;
            param_right = abs(self.positions(end) - self.positions(1)) * self.interface_right;
            
            % Calculate the deviation from the Kuprianov--Lukichev boundary
            % conditions, and store the results back into State instances
            s1.dg  = s1.dg  - ( eye(2) - s1.g*s0.gt )*N0*(  s0.g  - s1.g  )/param_left;
            s1.dgt = s1.dgt - ( eye(2) - s1.gt*s0.g )*Nt0*( s0.gt - s1.gt )/param_left;
            
            s2.dg  = s2.dg  - ( eye(2) - s2.g*s3.gt )*N3*(  s3.g  - s2.g  )/param_right;
            s2.dgt = s2.dgt - ( eye(2) - s2.gt*s3.g )*Nt3*( s3.gt - s2.gt )/param_right;

            % Vectorize the results of the calculations, and return it            
            residue = [s1.vectorize_dg s1.vectorize_dgt s2.vectorize_dg s2.vectorize_dgt]';
        end
        
        function gap = gap_calculate(self, position)
            % This function extracts the singlet components of the Green's
            % function at a given position, and then uses the gap equation
            % to calculate the superconducting gap at that point.
            
            singlets = zeros(size(self.energies));
            index    = self.position_index(position);
            
            % Extract the singlet components from the states
            for m=1:length(self.energies)
                singlets(m) = self.states(index,m).singlet;
            end
                
            % Create a cubic interpolation of the numerical data above,
            % multiplied by the tanh(ε/2T) kernel in the gap equation
            kernel = @(E) real(pchip(self.energies, singlets, E)) ...
                       .* tanh(E./(2*self.temperature));

            % Perform a numerical integration of the interpolation up to
            % the Debye cutoff (presumably the last element of 'energies')
            gap = self.scaling * integral(kernel, 0, self.energies(end));
        end
    end
end