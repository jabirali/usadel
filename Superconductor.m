% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-15
% Updated 2015-02-15


classdef Superconductor < handle
    % This defines a data structure that describes the physical state of
    % superconducting material for a range of positions and energies.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the internal variables for the data structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (GetAccess=public, SetAccess=public)
        positions   = [];
        energies    = [];
        gap         = [];
        states      = State.empty(0,0);
        
        bc_left     = State.empty(0);
        bc_right    = State.empty(0);
        temperature = 1e-16;
        scaling     = 1;
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
            self.states(length(positions), length(energies)) = State;
            for i=1:length(positions)
                for j=1:length(energies)
                    self.states(i,j) = Superconductor.Bulk(energies(j), 1);
                end
            end
            
            % Set the boundary conditions to empty states by default
            self.bc_left(length(energies))  = 0;
            self.bc_right(length(energies)) = 0;    
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
        
        function update_solution
        end
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods describing the behaviour of a superconductor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function dydx = jacobian(x,y)
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
                    - 2i*energy*g   ...
                    - gap*(SpinVector.Pauli.y - g * SpinVector.Pauli.y * g);
            
            d2gt =  - 2*dgt*N*g*dgt  ...
                    - 2i*energy*gt   ...
                    + gap*(SpinVector.Pauli.y - gt * SpinVector.Pauli.y * gt);
            
            % Fill the results of the calculations back into a 'State' object
            state.g   = dg;
            state.dg  = d2g;
            state.gt  = dgt;
            state.dgt = d2gt;
            
            % Pack the results into a state vector
            dydx = state.vectorize;
        end
        
        function boundary()
            % Calculates Kuprianov--Lukichev
        end
        
        
        
%         function gap = SuperconductingGap(states, energies, temperature, scaling)
%     % This function takes a the singlet component of the Green's function
%     % in a superconductor and the temperature of the system, and calculates
%     % the superconducting gap of the material.
%     %
%     % Input:   
%     %   states      This should be a vector containing either State objects
%     %               (numerical solution of the Usadel eq) or complex numbers
%     %               (the singlet components of the solutions).
%     %   energies    This should be a vector of energies which corresponds
%     %               to the states above.
%     %   temperature The temperature of the material.
%     %   cutoff      The Debye frequency cutoff of the material.
%     %   scaling     The superconducting gap is proportional to this
%     %               scaling constant (which is usually written N₀λ).
%     % 
%     % Output:  
%     %   gap         The superconducting gap (in general a complex number)
%     
%     % Extract the singlet components from the states, if necessary
%     singlets = zeros(1,numel(states));
%     if isa(states, 'State')
%         for n=1:length(singlets)
%             singlets(n) = states(n).singlet;
%         end
%     else
%         singlets = states;
%     end
% 
%     % Create a cubic interpolation of the numerical data above, multiplied
%     % by the tanh(ε/2T) kernel in the equation for the superconducting gap
%     kernel = @(E) pchip(energies, singlets, E) .* tanh(E./(2*temperature));
%     
%     % Perform a numerical integration of the interpolation up to the cutoff
%     gap    = scaling * integral(kernel, 0, energies(end));

    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define static methods (available without object instantiation)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    methods (Static)
        % Definition of static methods, which belong to the class and not the instance
        
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
    end
end