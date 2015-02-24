% This defines a data structure that describes the physical state of a
% ferromagnetic material with spin-orbit coupling for a given range
% of positions and energies. This class inherits the internal structure
% of the 'Metal' class.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Inspired by a similar program written by Sol Jacobsen
% Created 2015-02-16
% Updated 2015-02-19

classdef Ferromagnet < Metal
    properties (GetAccess=public, SetAccess=public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define the internal variables of the Ferromagnet class
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        exchange        = [0,0,0];           % Exchange field in the ferromagnet
        spinorbit       = SpinVector(0,0,0); % Spin-orbit field in the ferromagnet
    end
    
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods that instantiate the object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Ferromagnet(positions, energies, thouless, exchange, spinorbit)
            % This method constructs a Ferromagnet instance from a vector
            % of positions, a vector of energies, the Thouless energy, the
            % exchange field, and the spin-orbit field in the material.

            % Initialize the Metal superclass
            self@Metal(positions, energies, thouless);

            % Set the internal variables based on constructor arguments
            self.exchange  = exchange;
            self.spinorbit = spinorbit;
        end
        

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods which are useful for working with ferromagnets
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        
        function plot_comps(self)
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
        % Define methods that update the internal state of the object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                
        function update_coeff(self)
            % This function updates the vector of coefficients passed to
            % the functions 'jacobian' and 'boundary' when solving equations.
            
            % Coefficients in the equations for the Riccati parameter gamma
            self.coeff1{1} = -2;
            self.coeff1{2} = -2i/self.thouless;
            self.coeff1{3} = (self.spinorbit^2 - i*self.exchange*SpinVector.Pauli)/self.thouless;
            self.coeff1{4} = sqrt(2/self.thouless)    * self.spinorbit.x;
            self.coeff1{5} = sqrt(2/self.thouless)    * self.spinorbit.y;
            self.coeff1{6} = sqrt(2/self.thouless)    * self.spinorbit.z;
            self.coeff1{7} = (2i/sqrt(self.thouless)) * self.spinorbit.z;
            
            % Coefficients in the equations for the Riccati parameter gamma~
            self.coeff2{1} = self.coeff1{1};
            self.coeff2{2} = self.coeff1{2};
            self.coeff2{3} = conj(self.coeff1{3});
            self.coeff2{4} = conj(self.coeff1{4});
            self.coeff2{5} = conj(self.coeff1{5});
            self.coeff2{6} = conj(self.coeff1{6});
            self.coeff2{7} = conj(self.coeff1{7});
        end

        function update(self)
            % This function updates the internal state of the ferromagnet
            % by calling the other update methods. Always run this after
            % updating the boundary conditions, exchange field, or
            % spin-orbit field in the ferromagnet.
            
            % Update the state
            self.update_coeff;
            self.update_state;
    
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
            % This function takes a Ferromagnet object 'self', the
            % position 'x', the current state vector 'y', and an energy
            % as inputs, and calculates the Jacobian of the system. This 
            % is performed using the Riccati parametrized Usadel equations.
            
            % Extract the Riccati parameters and their derivatives
            [g,dg,gt,dgt] = State.unpack(y);
            
            % Calculate the normalization matrices
            N  = inv( eye(2) - g*gt );
            Nt = inv( eye(2) - gt*g );
            
            % Calculate the second derivatives of the Riccati parameters
            % according to the Usadel equations in the ferromagnet
            d2g  = self.coeff1{1} * dg*Nt*gt*dg                                                                  ...
                 + self.coeff1{2} * energy*g                                                                     ...
                 + (self.coeff1{3} * g - g * self.coeff2{3})                                                     ...
                 + (self.coeff1{4} * g + g * self.coeff2{4}) * Nt * (self.coeff2{4} + gt * self.coeff1{4} * g)   ...
                 + (self.coeff1{5} * g + g * self.coeff2{5}) * Nt * (self.coeff2{5} + gt * self.coeff1{5} * g)   ...
                 + (self.coeff1{6} * g + g * self.coeff2{6}) * Nt * (self.coeff2{6} + gt * self.coeff1{6} * g)   ...
                 + (self.coeff1{7} - g * self.coeff2{7} * gt) * N * dg                                           ...
                 + dg * Nt * (gt * self.coeff1{7} * g - self.coeff2{7});
                 
            d2gt = self.coeff2{1} * dgt*N*g*dgt                                                                  ...
                 + self.coeff2{2} * energy*gt                                                                    ...
                 + (self.coeff2{3} * gt - gt * self.coeff1{3})                                                   ...
                 + (self.coeff2{4} * gt + gt * self.coeff1{4}) * N * (self.coeff1{4} + g * self.coeff2{4} * gt)  ...
                 + (self.coeff2{5} * gt + gt * self.coeff1{5}) * N * (self.coeff1{5} + g * self.coeff2{5} * gt)  ...
                 + (self.coeff2{6} * gt + gt * self.coeff1{6}) * N * (self.coeff1{6} + g * self.coeff2{6} * gt)  ...
                 + (self.coeff2{7} - gt * self.coeff1{7} * g) * Nt * dgt                                         ...
                 + dgt * N * (g * self.coeff2{7} * gt - self.coeff1{7});
            
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
            dg1  = dg1  - ( eye(2) - g1*gt0 )*N0*(  g1  - g0  )/self.interface_left  + (self.coeff1{7} * g1  - g1  * self.coeff2{7})/2;
            dgt1 = dgt1 - ( eye(2) - gt1*g0 )*Nt0*( gt1 - gt0 )/self.interface_left  + (self.coeff2{7} * gt1 - gt1 * self.coeff1{7})/2;
            
            dg2  = dg2  - ( eye(2) - g2*gt3 )*N3*(  g2  - g3  )/self.interface_right + (self.coeff1{7} * g2  - g2  * self.coeff2{7})/2;
            dgt2 = dgt2 - ( eye(2) - gt2*g3 )*Nt3*( gt2 - gt3 )/self.interface_right + (self.coeff1{7} * g2  - g2  * self.coeff2{7})/2;
            
            % Vectorize the results of the calculations, and return it
            residue = State.pack(dg1,dgt1,dg2,dgt2);
        end        
    end
end
