% This defines a data structure that describes the physical state of a
% superconducting material for a given range of positions and energies.
% The most demanding calculations are parallellized using SPMD.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Inspired by a similar program written by Sol Jacobsen
% Created 2015-02-15
% Updated 2015-02-19

classdef Superconductor < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the internal variables for the data structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (GetAccess=public, SetAccess=public)
        % Properties that determine the physical characteristics of the system
        positions   = [];                    % Positions in the superconductor (relative to the material length)
        energies    = [];                    % Energies of the superconductor (relative to the Thouless energy)
        gap         = [];                    % Superconducting gap at each position
        states      = State.empty(0,0);      % Riccati parameters and their derivatives for each position and energy

        length          = 1;                 % Length of the system
        strength        = 1;                 % Material constant N₀λ
        temperature     = 0;                 % Temperature of the system (defaults to absolute zero)
        diffusion       = 1;                 % Diffusion constant
        interface_left  = inf;               % Interface parameter zeta (left)
        interface_right = inf;               % Interface parameter zeta (right)
        
        boundary_left   = State.empty(0);    % Boundary condition (left)   (default: vacuum)
        boundary_right  = State.empty(0);    % Boundary condition (right)  (default: vacuum)
        
        % Properties that determine the behavior of the program
        sim_error_abs = 1e-2;                % Maximum absolute error when simulating
        sim_error_rel = 1e-2;                % Maximum relative error when simulating
        sim_grid_size = 128;                 % Maximum grid size to use in simulations

        debug         = true;                % Whether to show intermediate results or not
        plot          = true;                % Whether to plot intermediate results or not
    end
    

    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define constructor and accessor methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

        function self = Superconductor(positions, energies, material_length, material_diffusion, material_strength)
            % Define a constructor which initializes the Superconductor
            % from a vector of positions, a vector of energies (where the
            % last element is the Debye cutoff), the length of the material
            % (scales the position vector), the material strength
            % (N₀λ in the gap equation), and the diffusion constant.
                        
            % Set default values based on constructor arguments
            self.positions = positions;
            self.energies  = energies;
            self.diffusion = material_diffusion;
            self.length    = material_length;
            self.strength  = material_strength;
            
            % Use eq. (3.34) in Tinkham to estimate the zero-temperature gap
            self.gap = abs(energies(end))/sinh(inv(material_strength)) ...
                     * ones(size(positions));
            
            % Set the maximum grid size for numerical calculations to 4x
            % the number of positions, rounded up to nearest power of two
            self.sim_grid_size = 2^(3+floor(log2(length(positions)-1)));
            
            % Initialize the internal state to a BCS bulk superconductor
            self.states(length(positions), length(energies)) = 0;
            for i=1:length(positions)
                for j=1:length(energies)
                    self.states(i,j) = Superconductor.Bulk(energies(j), self.gap(1));
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

            % Partially evaluate the Jacobian and boundary conditions
            % for the different superconductor energies, and put the
            % resulting lambda functions in a vector. These functions
            % are passed on to bvp6c when solving the equations.
            jc = {};
            bc = {};
            for m=1:length(self.energies)
                jc{m} = @(x,y) self.jacobian(self,x,y,self.energies(m));
                bc{m} = @(a,b) self.boundary(self,a,b,self.energies(m));
            end
                
            % Solving the differential equation is slow, and the solutions
            % at different energies are independent, so we parallelize this
            spmd
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
                    
                    % Try to solve the differential equation; if the solver
                    % returns an error, don't crash the program, but display a
                    % warning message and continue.
                    try
                        % Solve the differential equation, and evaluate the
                        % solution on the position vector of the superconductor 
                        solution = deval(bvp6c(jc{m},bc{m},initial,options), self.positions);

                        % Update the current state of the system based on the solution
                        for n=1:length(self.positions)
                            self.states(n,m) = State(solution(:,n));
                        end
                                    
                        % Progress information
                        self.print('[ %2.f / %2.f ]   iteration complete!', m, length(self.energies));

                    catch
                        % Display a warning message if the computation failed
                        self.print('[ %2.f / %2.f ]   iteration failed to converge, skipping...', m, length(self.energies));
                    end
                end
                
                % Small time delay to prevent the interpreter from getting sluggish
                pause(0.05);
            end
        end
        
        function update(self)
            % This function updates the internal state of the
            % superconductor object by calling first 'state_update' and 
            % then 'gap_update'. Always run this after changing the
            % boundary conditions or temperature of the system.
            
            % Update the state and gap estimates for the system
            gap_prev = mean(abs(self.gap));
            self.state_update;
            self.gap_update;
            gap_curr = mean(abs(self.gap));
            
            % If the relative change in gap was larger than 1%,
            % keep updating the estimate until we arrive at a
            % self-consistent solution of the problem
            while ~self.critical && (abs(1-gap_curr/gap_prev) > 1e-2) 
                % Status information
                self.print('Gap estimate changed by %2.f%%, recalculating...', 100*abs(gap_curr/gap_prev));

                % Update the estimate
                gap_prev = gap_curr;
                self.state_update;
                self.gap_update;
                gap_curr = mean(abs(self.gap));
                
                % Plot the current DOS
                if self.plot
                    self.plot_dos;
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define miscellaneous methods for showing the internal state
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


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define miscellaneous methods for use when debugging
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

        function print(self,varargin)
            % This function is used to print progress information if the
            % 'debug' flag is set to 'true'.
            
            if self.debug
                fprintf(':: Superconductor: %s\n', sprintf(varargin{:}));
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define static methods (available without object instantiation)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    methods (Static)        
        function result = Bulk(energy, gap)
            % This function takes as its input an energy and a superconducting gap,
            % and returns a State object with Riccati parametrized Green's functions
            % that corresponds to a BCS superconductor bulk state.
            theta = atanh(gap/(energy+0.001i));
            c     = cosh(theta);
            s     = sinh(theta);
            
            result = State([0,  s/(1+c); -s/(1+c), 0], 0, ...
                           [0, -s/(1+c);  s/(1+c), 0], 0);
        end
    
        function dydx = jacobian(self, x, y, energy)
            % This function takes a Superconductor object 'self', the
            % position 'x', the current state vector 'y', and an energy as
            % inputs, and calculates the Jacobian of the system. This is
            % performed using the Riccati parametrized Usadel equations.
            
            % Instantiate a 'State' object based on the state vector
            state = State(y);
            
            % Extract the superconducting gap
            gap = self.gap_lookup(x);
            
            % Calculate the Thouless energy
            thouless = self.diffusion/self.length^2;
            
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
                    - 2i*((energy+0.001i)/thouless)*g   ...
                    - (gap/thouless)*(SpinVector.Pauli.y - g * SpinVector.Pauli.y * g);
            
            d2gt =  - 2*dgt*N*g*dgt  ...
                    - 2i*((energy-0.001i)/thouless)*gt   ...
                    + (gap/thouless)*(SpinVector.Pauli.y - gt * SpinVector.Pauli.y * gt);
            
            % Fill the results of the calculations back into a 'State' object
            state.g   = dg;
            state.dg  = d2g;
            state.gt  = dgt;
            state.dgt = d2gt;
            
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
            
            % Calculate the deviation from the Kuprianov--Lukichev boundary
            % conditions, and store the results back into State instances
            s1.dg  = s1.dg  - ( eye(2) - s1.g*s0.gt )*N0*(  s0.g  - s1.g  )/self.interface_left;
            s1.dgt = s1.dgt - ( eye(2) - s1.gt*s0.g )*Nt0*( s0.gt - s1.gt )/self.interface_left;
            
            s2.dg  = s2.dg  - ( eye(2) - s2.g*s3.gt )*N3*(  s3.g  - s2.g  )/self.interface_right;
            s2.dgt = s2.dgt - ( eye(2) - s2.gt*s3.g )*Nt3*( s3.gt - s2.gt )/self.interface_right;

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
                   
            % Plot the current integral kernel
            if self.plot
                energies = linspace(0, self.energies(end), 100);
                plot(energies, kernel(energies));
                title('Integrand used in the gap equation');
                xlabel('Energy');
                ylabel('Re{f(E)} tanh(E/2T)');
                pause(0.05);
            end

            % Perform a numerical integration of the interpolation up to
            % the Debye cutoff (presumably the last element of 'energies')
            gap = -self.strength * integral(kernel, 0, self.energies(end));
        end
    end
end
