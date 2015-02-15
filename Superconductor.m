% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-15
% Updated 2015-02-15
%
% This defines a data structure that describes the physical state of
% superconducting material for a given range of positions and energies.


classdef Superconductor < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the internal variables for the data structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (GetAccess=public, SetAccess=public)
        positions   = [];
        energies    = [];
        gap         = [];
        states      = State.empty(0,0);
        
        temperature     = 1e-16;             % Temperature of the system
        scaling         = 1;                 % Material constant N₀λ
        diffusion       = 1;                 % Diffusion constant
        interface_left  = 1;                 % Interface parameter (left)
        interface_right = 1;                 % Interface parameter (right)
        
        boundary_left   = State.empty(0);    % Boundary condition (left)
        boundary_right  = State.empty(0);    % Boundary condition (right)
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the internal methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    methods
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
            
            % Set the boundary conditions to empty states by default
            self.boundary_left(length(energies))  = 0;
            self.boundary_right(length(energies)) = 0;    
        end
        
        function update_gap(self)
            % This function extracts the singlet component of the Green's
            % function of the superconductor at each position and energy,
            % and then uses the gap equation to update the current estimate
            % of the superconducting gap at equilibrium.
            
            singlets = zeros(size(self.energies));
            for n=1:length(self.positions)
                % Extract the singlet components from the states
                for m=1:length(self.energies)
                    singlets(m) = self.states(n,m).singlet;
                end
                
                % Create a cubic interpolation of the numerical data above,
                % multiplied by the tanh(ε/2T) kernel in the gap equation
                kernel = @(E) pchip(self.energies, singlets, E) ...
                           .* tanh(E./(2*self.temperature));

                % Perform a numerical integration of the interpolation up to
                % the Debye cutoff (presumably the last element of 'energies')
                self.gap(n) = self.scaling * integral(kernel, 0, self.energies(end));
            end
        end
        
        function update_state(self)
            % This function solves the Usadel equation numerically for the
            % given position and energy range, using the current stored 
            % estimate for the superconducting gap.

            % Set the accuracy of the numerical solution
            options = bvpset('AbsTol',1e-04,'RelTol',1e-04,'Nmax',1000);
            
            for m=1:length(self.energies)
                % Use the current state of the system as an initial guess
                current = self.states(:,m);
                for n=1:length(current)
                    current(n) = current(n).vectorize;
                end
                %@(n) self.states(n,m).vectorize;
                initial = bvpinit(self.positions, current);
                
                %
                jacobian = @(x,y) Superconductor.jacobian(x,y,self.energies(m));
                boundary = @(a,b) Superconductor.boundary(a,b);
                solution = bvp6c(jacobian,boundary,initial,options);
                
                % TODO:
                % Jacobian and Boundary: methods (require self) or static (need variables passing)
                % Implement the boundary conditions
                % 
            end
        end
        
        function update(self)
            % This function updates the internal state of the
            % superconductor by first solving the Usadel equation
            % numerically, and then calculating the superconducting gap
            % using the solution. Always run this after changing the
            % boundary conditions or temperature of the system.
            
            self.update_state;
            self.update_gap;
        end
        
        function result = critical(self)
            % This function returns whether or not the system is above
            % critical temperature, i.e. if the superconducting gap is zero
            % everywhere in the superconductor.
            result = ( max(abs(self.gap)) < 1e-4 );
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
    
        function dydx = jacobian(x, y, energy)
            % This function takes the position 'x' and current state vector 'y' as
            % inputs, and calculates the Jacobian of the system. This is performed
            % using the Riccati parametrized Usadel eqs in the *superconductor*.
            %
            % The function is nested, and can therefore access the variables of the
            % parent function to determine the energy and superconducting gap.
            
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
            % according to the Usadel equation in the superconductor
            d2g  =  - 2*dg*Nt*gt*dg ...
                    - 2i*(energy/diffusion)*g   ...
                    - (gap/diffusion)*(SpinVector.Pauli.y - g * SpinVector.Pauli.y * g);
            
            d2gt =  - 2*dgt*N*g*dgt  ...
                    - 2i*(energy/diffusion)*gt   ...
                    + (gap/diffusion)*(SpinVector.Pauli.y - gt * SpinVector.Pauli.y * gt);
            
            % Fill the results of the calculations back into a 'State' object
            state.g   = dg;
            state.dg  = d2g;
            state.gt  = dgt;
            state.dgt = d2gt;
            
            % Pack the results into a state vector
            dydx = state.vectorize;
        end
        
        function r = boundary(y1, y2)
            r = y1-y2;
        end
    end
end