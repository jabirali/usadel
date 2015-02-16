% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-16
% Updated 2015-02-16
%
% This defines a data structure that describes the physical state of
% ferromagnetic material with spin-orbit coupling for a given range
% of positions and energies.

classdef Ferromagnet < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the internal variables for the data structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (GetAccess=public, SetAccess=public)
        positions   = [];                    % Positions in the ferromagnet
        energies    = [];                    % Energies of the ferromagnet
        states      = State.empty(0,0);      % Green's functions for each position and energy

        exchange        = [0,0,0];           % Exchange field
        spinorbit       = SpinVector(0,0,0); % Spin-orbit field
        diffusion       = 1;                 % Diffusion constant
        interface_left  = inf;               % Interface parameter zeta (left)
        interface_right = inf;               % Interface parameter zeta (right)
        
        boundary_left   = State.empty(0);    % Boundary condition (left)   (default: vacuum)
        boundary_right  = State.empty(0);    % Boundary condition (right)  (default: vacuum)
                
        sim_error_abs = 1e-2;                % Maximum absolute error when simulating
        sim_error_rel = 1e-2;                % Maximum relative error when simulating
        sim_grid_size = 2048;                % Maximum grid size to use in simulations
    end
    

    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define constructor and accessor methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

        function self = Ferromagnet(positions, energies)
            % Define a constructor which initializes the Ferromagnet
            % from a vector of positions and a vector of energies
            self.positions = positions;
            self.energies  = energies;
            self.states(length(positions), length(energies)) = 0;
            
            % Initialize the internal state to a weak bulk superconductor
            self.states(length(positions), length(energies)) = 0;
            for i=1:length(positions)
                for j=1:length(energies)
                    self.states(i,j) = Superconductor.Bulk(energies(j), 0.1);
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
            
%            tic;
%            disp(sprintf('\n:: Ferromagnet: updating state...'));
            for m=1:length(self.energies)
                % Vectorize the current state of the system for the given
                % energy, and use it as an initial guess for the solution
                current = zeros(16,length(self.positions));
                for n=1:length(self.positions)
                    current(:,n) = self.states(n,m).vectorize;
                end
                initial = bvpinit(self.positions', current);
                
                % Partially evaluate the Jacobian and boundary conditions
                % for the current energy
                jc = @(x,y) Ferromagnet.jacobian(self,x,y,self.energies(m));
                bc = @(a,b) Ferromagnet.boundary(self,a,b,self.energies(m));
                
                % Solve the differential equation, and evaluate the
                % solution on the position vector of the ferromagnet 
                solution = deval(bvp6c(jc,bc,initial,options), self.positions);

                % Update the current state of the system based on the solution
                for n=1:length(self.positions)
                    self.states(n,m) = State(solution(:,n));
                end
                
                % Progress information
 %               disp(sprintf('--     ETA: %2.f sec.', toc*(1-m/length(self.energies))));
            end
        end
        
        function update(self)
            % This function updates the internal state of the ferromagnet
            % by calling the 'update_*' methods. Always run this after
            % changing the boundary conditions or parameters of the system.
            
            self.state_update;
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
            
            % Extract diffusion constant and background fields
            diffusion = self.diffusion;
            exchange  = self.exchange;
            spinorbit = self.spinorbit;
            
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
            d2g  =  - 2*dg*Nt*gt*dg                                                             ...
                    - 2i*((energy+0.001i)/diffusion)*g                                                   ...
                    - i*(exchange/diffusion)*(SpinVector.Pauli*g - g*conj(SpinVector.Pauli))    ...
                    + 2i*(spinorbit.z + g*conj(spinorbit.z)*gt)*N*dg                            ...
                    + 2i*dg*Nt*(conj(spinorbit.z) + gt*spinorbit.z*g)                           ...
                    + 2*(spinorbit*g + g*conj(spinorbit))*Nt*(conj(spinorbit) + gt*spinorbit*g) ...
                    + (spinorbit^2*g - g*conj(spinorbit)^2);
            
            d2gt =  - 2*dgt*N*g*dgt                                                             ...
                    - 2i*((energy+0.001i)/diffusion)*gt                                                  ...
                    + i*(exchange/diffusion)*(conj(SpinVector.Pauli)*gt - gt*SpinVector.Pauli)  ...
                    - 2i*(conj(spinorbit.z) + gt*spinorbit.z*g)*Nt*dgt                          ...
                    - 2i*dgt*N*(spinorbit.z + g*conj(spinorbit.z)*gt)                           ...
                    + 2*(conj(spinorbit)*gt + gt*spinorbit)*N*(spinorbit + g*conj(spinorbit)*gt)...
                    + (conj(spinorbit)^2*gt - gt*spinorbit^2);
            
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
             
            % Calculate the normalization matrices
            N0  = inv( eye(2) - s0.g*s0.gt );
            Nt0 = inv( eye(2) - s0.gt*s0.g );

            N1  = inv( eye(2) - s1.g*s1.gt );
            Nt1 = inv( eye(2) - s1.gt*s1.g );

            N2  = inv( eye(2) - s2.g*s2.gt );
            Nt2 = inv( eye(2) - s2.gt*s2.g );

            N3  = inv( eye(2) - s3.g*s3.gt );
            Nt3 = inv( eye(2) - s3.gt*s3.g );
            
            % Calculate interface parameters in the Kuprianov-Lukichev B.C.
            param_left  = abs(self.positions(end) - self.positions(1)) * self.interface_left;
            param_right = abs(self.positions(end) - self.positions(1)) * self.interface_right;
            
            % Calculate the deviation from the Kuprianov--Lukichev boundary
            % conditions, and store the results back into State instances
            s1.dg  = s1.dg  - ( eye(2) - s1.g*s0.gt )*N0*(  s1.g  - s0.g  )/param_left  ...
                   - 1i * (self.spinorbit.z * s1.g + s1.g * conj(self.spinorbit.z));
            s1.dgt = s1.dgt - ( eye(2) - s1.gt*s0.g )*Nt0*( s1.gt - s0.gt )/param_left  ...
                   + 1i * (conj(self.spinorbit.z) * s1.gt + s1.gt * self.spinorbit.z);
            
            s2.dg  = s2.dg  - ( eye(2) - s2.g*s3.gt )*N3*(  s2.g  - s3.g  )/param_right ...
                   - 1i * (self.spinorbit.z * s2.g + s2.g * conj(self.spinorbit.z));
            s2.dgt = s2.dgt - ( eye(2) - s2.gt*s3.g )*Nt3*( s2.gt - s3.gt )/param_right ...
                   + 1i * (conj(self.spinorbit.z) * s2.gt + s2.gt * self.spinorbit.z);

            % Vectorize the results of the calculations, and return it            
            residue = [s1.vectorize_dg s1.vectorize_dgt s2.vectorize_dg s2.vectorize_dgt]';
        end
    end
end