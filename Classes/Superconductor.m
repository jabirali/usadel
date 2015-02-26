% This defines a data structure that describes the physical state of a
% superconducting material for a given range of positions and energies.
% This class inherits the internal structure of the 'Metal' class.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Inspired by a similar program written by Sol Jacobsen
% Created 2015-02-15
% Updated 2015-02-26

classdef Superconductor < Metal
    properties (GetAccess=public, SetAccess=public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define the internal variables of the Superconductor class
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        temperature = 0;                                % Temperature of the system (relative to the critical temperature of a bulk superconductor)
        strength    = 1;                                % Strength of the superconductor (the material constant N0V which appears in the gap equation)
        gap         = griddedInterpolant([0,1],[1,1]);  % Superconducting gap as a function of position (relative to the zero-temperature gap of a bulk superconductor)
        
    end
    
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods that instantiate the object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Superconductor(positions, energies, thouless, strength)
            % This method constructs a Superconductor instance from a vector
            % of positions, a vector of energies, the Thouless energy, and
            % the strength of the superconductivity (material constant N0V).

            % Initialize the Metal superclass
            self@Metal(positions, energies, thouless);

            % Set the internal variables based on constructor arguments
            self.strength = strength;
        end

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods which are useful for working with superconductors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function result = critical(self)
            % This function returns whether or not the system is above
            % critical temperature. The superconductor is considered
            % critical if the maximal superconducting gap is below a
            % certain threshould value.
            
            result = ( abs(self.mean_gap) < 1e-4 );
        end        
        
        function result = mean_gap(self)
            % This function returns the mean superconducting gap.
            
            result = mean(self.gap.Values);
        end
        
        function plot_gap(self)
            % Plot the superconducting gap as a function of position
            xs = linspace(self.positions(1), self.positions(end));
            plot(xs, self.gap(xs));
            xlabel('Relative position');
            ylabel('Superconducting gap');
        end        

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods that update the internal state of the object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        
        function update_gap(self)
            % This function updates the gap function, which contains the 
            % current estimate of the superconducting gap in the material.
            
            % Calculate the gap at the positions in the system
            gaps = [];
            for n=1:length(self.positions)
                gaps(n) = self.calculate_gap(self, self.positions(n));
            end
            
            % Create a piecewise cubic interpolation of the results
            self.gap = griddedInterpolant(self.positions, gaps, 'pchip');
        end
        
        function update_coeff(self)
            % This function updates the vector of coefficients passed to
            % the functions 'jacobian' and 'boundary' when solving equations.
            
            % Coefficients in the equations for the Riccati parameter gamma
            self.coeff1{1} = -2i/self.thouless;
            self.coeff1{2} = SpinVector.Pauli.y/self.thouless;
            
            % Coefficients in the equations for the Riccati parameter gamma~
            self.coeff2{1} =  self.coeff1{1};
            self.coeff2{2} = -self.coeff1{2};
        end

        function update(self)
            % This function updates the internal state of the superconductor
            % by calling the other update methods. Always run this after
            % updating the boundary conditions or physical properties of
            % the superconductor.
            
            % Update the state
            self.update_coeff;
            self.update_state;
            self.update_gap;
                
            % Plot the current DOS
            if self.plot
                self.plot_dos;
            end
        end
    end        
    
    methods (Static)        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define static methods (available without object instantiation)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

        function dydx = jacobian(self, x, y, energy)
            % This function takes a Superconductor object 'self', the
            % position 'x', the current state vector 'y', and an energy
            % as inputs, and calculates the Jacobian of the system. This 
            % is performed using the Riccati parametrized Usadel equations.
            
            % Extract the Riccati parameters and their derivatives
            [g,dg,gt,dgt] = State.unpack(y);
            
            % Retrieve the superconducting gap at this point
            gap = self.gap(x);
            
            % Calculate the normalization matrices
            N  = inv( eye(2) - g*gt );
            Nt = inv( eye(2) - gt*g );
            
            % Calculate the second derivatives of the Riccati parameters
            % according to the Usadel equations in the superconductor
            d2g  = -2 * dg*Nt*gt*dg + self.coeff1{1} * energy*g         ...
                 + gap * (self.coeff1{2} + g*self.coeff2{2}*g);
             
            d2gt = -2 * dgt*N*g*dgt + self.coeff2{1} * energy*gt        ...
                 + gap * (self.coeff2{2} + gt*self.coeff1{2}*gt);
            
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
            dg1  = dg1  - ( eye(2) - g1*gt0 )*N0*(  g0  - g1  )/self.interface_left;
            dgt1 = dgt1 - ( eye(2) - gt1*g0 )*Nt0*( gt0 - gt1 )/self.interface_left;
            
            dg2  = dg2  - ( eye(2) - g2*gt3 )*N3*(  g3  - g2  )/self.interface_right;
            dgt2 = dgt2 - ( eye(2) - gt2*g3 )*Nt3*( gt3 - gt2 )/self.interface_right;
            
            % Vectorize the results of the calculations, and return it
            residue = State.pack(dg1,dgt1,dg2,dgt2);
        end
        
        function gap = calculate_gap(self, position)
            % This function extracts the singlet components of the Green's
            % function at a given position, and then uses the gap equation
            % to calculate the superconducting gap at that point.
            
            % Extract the singlet components from the states
            singlets = zeros(size(self.energies));
            index    = self.position_index(position);
            for m=1:length(self.energies)
                singlets(m) = self.states(index,m).singlet;
            end
            
            % Create a cubic interpolation of the numerical data above,
            % multiplied by the tanh(E/2T) kernel in the gap equation.
            % Using eq. (4.34) in Fossheim & Sudbø, we can rewrite the 
            % argument of the hyperbolic tangent such that E is measured 
            % relative to the zero-temperature gap, and T relative to Tc.
            kernel = @(E) real(pchip(self.energies, singlets, E)) ...
                       .* tanh(0.881939 * E/self.temperature);
                   
            % Perform a numerical integration of the interpolation up to
            % the Debye cutoff. The Debye cutoff has been calculated from
            % the superconductor strength using eq. (3.34) in Tinkham.
            gap = -self.strength * integral(kernel, 0, sinh(inv(self.strength)));
        end
        
        function result = Bulk(energy, gap)
            % This function takes as its input an energy and a superconducting gap,
            % and returns a State object with Riccati parametrized Green's functions
            % that corresponds to a BCS superconductor bulk state. We also
            % include an infinitesimal triplet contribution.
            theta = atanh(gap/(energy+1e-16i));
            c     = cosh(theta);
            s     = sinh(theta);
            
            result = State([0,  s/(1+c); -s/(1+c), 0], 0, ...
                           [0, -s/(1+c);  s/(1+c), 0], 0);
        end
    end
end
