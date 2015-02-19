% This defines a data structure that describes the physical state of
% ferromagnetic material with spin-orbit coupling for a given range
% of positions and energies. The most demanding calculations are
% parallellized using SPMD.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Inspired by a similar program written by Sol Jacobsen
% Created 2015-02-16
% Updated 2015-02-19

classdef Ferromagnet < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the internal variables for the data structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (GetAccess=public, SetAccess=public)
        % Properties that determine the physical characteristics of the system
        positions   = [];                    % Positions in the ferromagnet (relative to the material length)
        energies    = [];                    % Energies of the ferromagnet (relative to the Thouless energy)
        states      = State.empty(0,0);      % Green's functions for each position and energy

        length          = 1;                 % Length of the system
        exchange        = [0,0,0];           % Exchange field
        spinorbit       = SpinVector(0,0,0); % Spin-orbit field
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
        plot          = true;                % Whether to plot inter
    end
    

    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define constructor and accessor methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

        function self = Ferromagnet(positions, energies, material_length, material_diffusion)
            % Define a constructor which initializes the Ferromagnet
            % from a vector of positions, a vector of energies (where the
            % last element is the Debye cutoff), the length of the material
            % (scales the position vector), and the diffusion constant.
            
            % Set default values based on constructor arguments
            self.positions = positions;
            self.energies  = energies;
            self.diffusion = material_diffusion;
            self.length    = material_length;
            self.states(length(positions), length(energies)) = 0;

            % Set the maximum grid size for numerical calculations to 4x
            % the number of positions, rounded up to nearest power of two
            self.sim_grid_size = 2^(3+floor(log2(length(positions)-1)));

            % Initialize the internal state to a bulk superconductor
            self.states(length(positions), length(energies)) = 0;
            for i=1:length(positions)
                for j=1:length(energies)
                    self.states(i,j) = Superconductor.Bulk(energies(j), 1);
                end
            end
            
            % Set the boundary conditions to empty states by default
            self.boundary_left(length(energies))  = 0;
            self.boundary_right(length(energies)) = 0;    
        end
        
        function index = position_index(self, position)
            % Returns the vector index corresponding to a given energy value
            index = find(abs(self.positions-position) < 1e-5, 1, 'first');
            if ~isscalar(index)
                error('Ferromagnet.position_index: Provided value is not in the position vector!');
            end
        end

        function index = energy_index(self, energy)
            % Returns the vector index corresponding to a given energy value
            index = find(abs(self.energies-energy) < 1e-5, 1, 'first');
            if ~isscalar(index)
                error('Ferromagnet.energy_index: Provided value is not in the energy vector!');
            end
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods that update the internal state of the object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        
        function state_update(self)
            % This function solves the Usadel equation numerically for the
            % given position and energy range, using the current stored 
            % exchange field and spin-orbit field.
            
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
                        % solution on the position vector of the ferromagnet 
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
            % This function updates the internal state of the ferromagnet
            % by calling the 'update_*' methods. Always run this after
            % changing the boundary conditions or parameters of the system.
            
            self.state_update;
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
            srtc    = zeros(length(self.positions), 1);
            lrtc    = zeros(length(self.positions), 1);
            for m=1:length(self.energies)
                for n=1:length(self.positions)
                    singlet(n) = singlet(n) + norm(self.states(n,m).singlet);
                    triplet(n) = triplet(n) + norm(self.states(n,m).triplet);
                    srtc(n) = srtc(n) + norm(self.states(n,m).srtc(self.exchange));
                    lrtc(n) = lrtc(n) + norm(self.states(n,m).lrtc(self.exchange));
                end
            end
                        
            % Plot cubic interpolations of the results
            positions = linspace(0, self.positions(end), 100);
            if norm(self.exchange) == 0
                % If there is no exchange field, don't distinguish SRT/LRT
                plot(positions, pchip(self.positions, singlet, positions),  ...
                     positions, pchip(self.positions, triplet, positions));
                legend('Singlet', 'Triplet');
            else
                % If there is an exchange field, distinguish SRT/LRT
                plot(positions, pchip(self.positions, singlet, positions),  ...
                     positions, pchip(self.positions, srtc,    positions),  ...
                     positions, pchip(self.positions, lrtc,    positions));
                legend('Singlet', 'Short-Range Triplet', 'Long-Range Triplet');
            end
            xlabel('Relative position');
            ylabel('Distribution');
        end
 
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define miscellaneous methods for use when debugging
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

        function print(self,varargin)
            % This function is used to print progress information if the
            % 'debug' flag is set to 'true'.
            
            if self.debug
                fprintf(':: Ferromagnet:    %s\n', sprintf(varargin{:}));
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define static methods (available without object instantiation)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    methods (Static)            
        function dydx = jacobian(self, x, y, energy)
            % This function takes a Ferromagnet object 'self', the
            % position 'x', the current state vector 'y', and an energy as
            % inputs, and calculates the Jacobian of the system. This is
            % performed using the Riccati parametrized Usadel equations.
            
            % Instantiate a 'State' object based on the state vector
            state = State(y);
                        
            % Introduce some short notations
            h = self.exchange;
            A = self.spinorbit;
            L = self.length;
            D = self.diffusion;

            % Extract the Riccati parameters and their derivatives
            g   = state.g;
            dg  = state.dg;
            gt  = state.gt;
            dgt = state.dgt;
            
            % Calculate the normalization matrices
            N  = inv( eye(2) - g*gt );
            Nt = inv( eye(2) - gt*g );
            
            % Calculate the second derivatives of the Riccati parameters
            % according to the Usadel equation in the ferromagnet
            d2g  =  - 2*dg*Nt*gt*dg                                                       ...
                    - (2i*L^2/D) * (energy+0.001i) * g                                    ...
                    - (1i*L^2/D) * h * (SpinVector.Pauli*g - g*conj(SpinVector.Pauli))    ...
                    + (2i*L/D)   * (A.z + g*conj(A.z)*gt)*N*dg                            ...
                    + (2i*L/D)   * dg*Nt*(conj(A.z) + gt*A.z*g)                           ...
                    + (2*L^2/D)  * (A*g + g*conj(A))*Nt*(conj(A) + gt*A*g)                ...
                    + (L^2/D)    * (A^2*g - g*conj(A)^2);
            
            d2gt =  - 2*dgt*N*g*dgt                                                       ...
                    - (2i*L^2/D) * (energy-0.001i) * gt                                   ...
                    + (1i*L^2/D) * h * (conj(SpinVector.Pauli)*gt - gt*SpinVector.Pauli)  ...
                    - (2i*L/D)   * (conj(A.z) + gt*A.z*g)*Nt*dgt                          ...
                    - (2i*L/D)   * dgt*N*(A.z + g*conj(A.z)*gt)                           ...
                    + (2*L^2/D)  * (conj(A)*gt + gt*A)*N*(A + g*conj(A)*gt)               ...
                    + (L^2/D)    * (conj(A)^2*gt - gt*A^2);
            
            % Fill the results of the calculations back into a 'State' object
            state.g   = dg;
            state.dg  = d2g;
            state.gt  = dgt;
            state.dgt = d2gt;
            
            % Pack the results into a state vector
            dydx = state.vectorize;
        end
        
        function residue = boundary(self, y1, y2, energy)
            % This function takes a Ferromagnet object 'self', the
            % position 'x', the current state vector 'y', and an energy as
            % inputs, and calculates the Kuprianov-Lukichev boundary
            % conditions for the system. 

            % State in the material to the left of the ferromagnet
            s0   = State(self.boundary_left(self.energy_index(energy)));
            
            % State at the left end of the ferromagnet
            s1   = State(y1);
            
            % State at the right end of the ferromagnet
            s2   = State(y2);
            
            % State in the material to the right of the ferromagnet
            s3   = State(self.boundary_right(self.energy_index(energy)));
            
            % The length of the system
            L = self.length;
             
            % Calculate the normalization matrices
            N0  = inv( eye(2) - s0.g*s0.gt );
            Nt0 = inv( eye(2) - s0.gt*s0.g );

            N1  = inv( eye(2) - s1.g*s1.gt );
            Nt1 = inv( eye(2) - s1.gt*s1.g );

            N2  = inv( eye(2) - s2.g*s2.gt );
            Nt2 = inv( eye(2) - s2.gt*s2.g );

            N3  = inv( eye(2) - s3.g*s3.gt );
            Nt3 = inv( eye(2) - s3.gt*s3.g );
                        
            % Calculate the deviation from the Kuprianov--Lukichev boundary
            % conditions, and store the results back into State instances
            s1.dg  = s1.dg  - ( eye(2) - s1.g*s0.gt )*N0*(  s1.g  - s0.g  )/self.interface_left  ...
                   - (1i*L) * (self.spinorbit.z * s1.g + s1.g * conj(self.spinorbit.z));
            s1.dgt = s1.dgt - ( eye(2) - s1.gt*s0.g )*Nt0*( s1.gt - s0.gt )/self.interface_left  ...
                   + (1i*L) * (conj(self.spinorbit.z) * s1.gt + s1.gt * self.spinorbit.z);
            
            s2.dg  = s2.dg  - ( eye(2) - s2.g*s3.gt )*N3*(  s2.g  - s3.g  )/self.interface_right ...
                   - (1i*L) * (self.spinorbit.z * s2.g + s2.g * conj(self.spinorbit.z));
            s2.dgt = s2.dgt - ( eye(2) - s2.gt*s3.g )*Nt3*( s2.gt - s3.gt )/self.interface_right ...
                   + (1i*L) * (conj(self.spinorbit.z) * s2.gt + s2.gt * self.spinorbit.z);

            % Vectorize the results of the calculations, and return it            
            residue = [s1.vectorize_dg s1.vectorize_dgt s2.vectorize_dg s2.vectorize_dgt]';
        end
    end
end