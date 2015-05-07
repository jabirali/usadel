% This defines a data structure that describes the physical state of a
% superconducting material for a given range of positions and energies.
% This class inherits the internal structure of the 'Metal' class.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Inspired by a similar program written by Sol Jacobsen
% Created 2015-02-15
% Updated 2015-05-05

classdef Superconductor < Metal
    properties (GetAccess=public, SetAccess=public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define the internal variables of the Superconductor class
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        temperature = 0;                                % Temperature of the system (relative to the critical temperature of a bulk superconductor)
        strength    = 0.2;                              % Strength of the superconductor (the material constant N0V which appears in the gap equation)
        gap         = griddedInterpolant([0,1],[1,1]);  % Superconducting gap as a function of position (relative to the zero-temperature gap of a bulk superconductor)
        phase       = griddedInterpolant([0,1],[0,0]);  % Superconducting phase as a function of position (defaults to zero)
        complex     = false;                            % Whether we use a gauge where the superconducting phase is zero
        locked      = false;                            % Whether the superconducting gap and phase are locked to constant values or not
    end
    
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods that instantiate the object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Superconductor(positions, energies, thouless, strength)
            % This method constructs a Superconductor instance from a vector
            % of positions, a vector of energies, the Thouless energy, and
            % the strength of the superconductivity (material constant N0V).
            %
            % Note: for self-consistent solutions to work properly, the energy
            %       vector should extend up to the Debye cutoff cosh(1/N0V).

            % Initialize the Metal superclass
            self@Metal(positions, energies, thouless);

            % Set the internal variables based on constructor arguments
            self.strength = strength;
        end

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods which are useful for working with superconductors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function gap_reduce(self)
            % This function reduces the current value of the superconducting
            % gap everywhere in the material. (This can be useful for
            % accelerating the convergence of programs that search
            % for the critical temperature of a hybrid system.)
            
            self.gap.Values = self.gap.Values./2;
            self.update_coeff;
            self.update_state;
            self.update_gap;
        end        
        
        function result = gap_left(self)
            % This function returns the left superconducting gap.
            
            result = abs(self.gap.Values(1));
        end

        function result = gap_right(self)
            % This function returns the right superconducting gap.
            
            result = abs(self.gap.Values(end));
        end

        function result = gap_mean(self)
            % This function returns the mean superconducting gap.
            
            result = mean(abs(self.gap.Values));
        end
        
        function result = gap_max(self)
            % This function returns the maximal superconducting gap.
            
            result = max(abs(self.gap.Values));
        end
        
        function result = gap_min(self)
            % This function returns the minimal superconducting gap.
            
            result = min(abs(self.gap.Values));
        end
        
        function phase_set(self, phase)
            % This function sets the superconducting phase at all positions,
            % and updates the internal state of the superconductor to a bulk
            % material with the correct phase.
            
            % Set the superconducting phase
            self.phase = griddedInterpolant([0,1],[phase,phase]);
            
            % Set the internal state to a bulk superconductor with a phase
            for i=1:length(self.positions)
                for j=1:length(self.energies)
                    self.states(i,j) = Superconductor.Bulk(self.energies(j), 1, phase);
                end
            end
        end
        
        function plot_gap(self)
            % Plot the superconducting gap as a function of position
            xs = linspace(self.positions(1), self.positions(end));
            plot(xs, self.gap(xs));
            xlabel('Relative position');
            ylabel('Superconducting gap');
        end        

        function plot_phase(self)
            % Plot the superconducting phase as a function of position
            xs = linspace(self.positions(1), self.positions(end));
            plot(xs, self.phase(xs)/pi);

            axis([0 1 -1 1]);
            set(gca, 'XTick', [0,1/4,1/2,3/4,1]);
            set(gca, 'YTick', [-1,-1/2,0,1/2,1]);

            xlabel('Relative position');
            ylabel('Superconducting phase (\pi)');
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods that update the internal state of the object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        
        function update_gap(self)
            % This function updates the gap function, which contains the 
            % current estimate of the superconducting gap in the material.
            
            if ~self.locked
                % Calculate the gap and phase at every position in the system
                gaps   = zeros(1,length(self.positions));
                phases = zeros(1,length(self.positions));
                for n=1:length(self.positions)
                    [gaps(n),phases(n)] = self.calculate_gap(self, self.positions(n));
                end
            
                % Create a piecewise cubic interpolation of the results
                self.gap   = griddedInterpolant(self.positions, gaps,   'pchip');
                self.phase = griddedInterpolant(self.positions, phases, 'pchip');
            end
        end
        
        function update_coeff(self)
            % This function updates the vector of coefficients passed to
            % the functions 'jacobian' and 'boundary' when solving equations.
            
            % Call the standard 'Metal' version of the method
            update_coeff@Metal(self);
            
            % Coefficients in the equations for the Riccati parameter gamma
            self.coeff1{1} = -2i/self.thouless;
            self.coeff1{2} = -SpinVector.Pauli.y/self.thouless;
            
            % Coefficients in the equations for the Riccati parameter gamma~
            self.coeff2{1} =  self.coeff1{1};
            self.coeff2{2} = -self.coeff1{2};
        end

        function update(self)
            % This function updates the internal state of the superconductor
            % by calling the other update methods. Always run this after
            % updating the boundary conditions or physical properties of
            % the superconductor, or after restoring from a backup.
            
            % Update the state
            self.update_gap;
            self.update_coeff;
            self.update_state;
            self.update_gap;
                
            % Plot the current DOS
            if self.plot
                self.plot_dos_surf;
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
            gap  = self.gap(x) * exp(i*self.phase(x));
            gapt = conj(gap);
            
            % Calculate the normalization matrices
            N  = inv( eye(2) - g*gt );
            Nt = inv( eye(2) - gt*g );
            
            % Calculate the second derivatives of the Riccati parameters
            % according to the Usadel equations in the superconductor
            d2g  = -2 * dg*Nt*gt*dg + self.coeff1{1} * (energy+1e-3i)*g ...
                 + gap * self.coeff1{2} + gapt * g*self.coeff2{2}*g;
             
            d2gt = -2 * dgt*N*g*dgt + self.coeff2{1} * (energy+1e-3i)*gt...
                 + gapt * self.coeff2{2} + gap * gt*self.coeff1{2}*gt;
            
            % Pack the results into a state vector
            dydx = State.pack(dg,d2g,dgt,d2gt);
        end
        
        function [gap,phase] = calculate_gap(self, position)
            % This function returns a vector [gap,phase] with the
            % superconducting gap and phase at the given position.
            % This is done by invoking either 'calculate_gap_real'
            % or 'calculate_gap_complex', depending on whether 
            % the parameter 'complex' is true or false.
            
            if self.complex
                    [gap,phase] = self.calculate_gap_complex(self, position);
            else
                    [gap,phase] = self.calculate_gap_real(self, position);
            end
        end
        
        function [result_gap,result_phase] = calculate_gap_real(self, position)
            % This function extracts the singlet components of the Green's
            % function at a given position, and then uses the gap equation
            % to calculate the superconducting gap at that point.
            % [This version works for purely real values of the gap.]
            
            % Extract the singlet components from the states
            singlets = zeros(size(self.energies));
            index    = self.position_index(position);
            for m=1:length(self.energies)
                singlets(m) = real(self.states(index,m).singlet);
            end
            
            % Create a cubic interpolation of the numerical data above
            singlet = griddedInterpolant(self.energies, singlets, 'pchip');
            
            % This is the expression for the gap equation integrand. Using
            % eq. (4.34) in Fossheim & Sudbø, we have rewritten the argument
            % of tanh(E/2T) such that E is measured relative to the
            % zero-temperature gap, while T is measured relative to Tc.
            kernel = @(E) singlet(E) .* tanh(0.881939 * E/self.temperature);
                   
            % Perform a numerical integration of the interpolation up to
            % the Debye cutoff. The Debye cutoff has been calculated from
            % the superconductor strength using eq. (3.34) in Tinkham.
            % The reason for cosh(1/N0V) instead of sinh(1/N0V), is that we
            % integrate the quasiparticle energy and not the kinetic energy.
            result_gap   = self.strength * integral(kernel, self.energies(1), cosh(1/self.strength));
            result_phase = 0;
        end
        
        
        function [result_gap,result_phase] = calculate_gap_complex(self, position)
            % This function uses the gap equation to calculate the super-
            % conducting gap and superconducting phase at a given position.
            % [This version works for general complex values of the gap.]
            
            % Calculate the singlet component of the anomalous Green's
            % functions and its tilde conjugate, and then their difference.
            index    = self.position_index(position);
            singlets = zeros(size(self.energies));
            for m=1:length(self.energies)
                singlets(m) = (self.states(index,m).singlet - conj(self.states(index,m).singlett))/2;
            end

            % Create cubic interpolations of the numerical data above
            singletR = griddedInterpolant(self.energies, real(singlets), 'pchip');
            singletI = griddedInterpolant(self.energies, imag(singlets), 'pchip');
            
            % This is the expression for the gap equation integrand. Using
            % eq. (4.34) in Fossheim & Sudbø, we have rewritten the argument
            % of tanh(E/2T) such that E is measured relative to the
            % zero-temperature gap, while T is measured relative to Tc.
            kernelR = @(E) singletR(E) .* tanh(0.881939 * E/self.temperature);
            kernelI = @(E) singletI(E) .* tanh(0.881939 * E/self.temperature);
                   
            % Perform a numerical integration of the interpolation up to
            % the Debye cutoff. The Debye cutoff has been calculated from
            % the superconductor strength using eq. (3.34) in Tinkham.
            % The reason for cosh(1/N0V) instead of sinh(1/N0V), is that we
            % integrate the quasiparticle energy and not the kinetic energy.
            gapR = self.strength * integral(kernelR, self.energies(1), cosh(1/self.strength));
            gapI = self.strength * integral(kernelI, self.energies(1), cosh(1/self.strength));
         
            % Extract the superconducting gap and phase from the results
            result       = gapR + 1i*gapI;
            result_gap   = norm(result);
            result_phase = phase(result);
        end
        
        function result = Bulk(energy, gap, phase)
            % This function takes as its input an energy and a superconducting
            % gap, and returns a State object with Riccati parametrized Green's
            % functions that correspond to a BCS superconductor bulk state.
                        
            theta = atanh(gap/(energy+1e-3i));
            c     = cosh(theta);
            s     = sinh(theta);
            g     = s/(1+c);
            
            result = State([0,  g; -g, 0]*exp( 1i*phase), 0, ...
                           [0, -g;  g, 0]*exp(-1i*phase), 0);
        end
    end
end
