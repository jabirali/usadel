% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-15
% Updated 2015-02-15
%
% This defines a data structure that describes the physical state of
% ferromagnetic material with spin-orbit coupling for a given range
% of positions and energies.


classdef Ferromagnet < handle   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the internal variables for the data structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (GetAccess=public, SetAccess=public)
        positions   = [];
        energies    = [];
        exchange    = [0,0,0];
        spinorbit   = SpinVector(0,0,0);
        states      = State.empty(0,0);
        
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
        function self = Ferromagnet(positions, energies)
            % Define a constructor which initializes the Ferromagnet
            % from a vector of positions and a vector of energies
            self.positions = positions;
            self.energies  = energies;
            self.states(length(positions), length(energies)) = 0;
            
            % Set the boundary conditions to empty states by default
            self.boundary_left(length(energies))  = 0;
            self.boundary_right(length(energies)) = 0;    
        end
        
        function update_state(self)
            % This function solves the Usadel equation numerically in the
            % ferromagnet, using the stored exchange and spin-orbit fields.
            
            % TODO
        end
        
        function update(self)
            % This function updates the internal state of the ferromagnet
            % by solving the Usadel equation numerically Always run this
            % after changing the boundary conditions, exchange field, or
            % spin-orbit field of the system.
            
            self.update_state;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define static methods (available without object instantiation)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    methods (Static)         
        function dydx = jacobian(x,y)
            % This function takes the position 'x' and current state vector 'y' as
            % inputs, and calculates the Jacobian of the system. This is performed
            % using the Riccati parametrized Usadel eqs in the *ferromagnet*.
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
            d2g  =  - 2*dg*Nt*gt*dg                                                             ...
                    - 2i*(energy/diffusion)*g                                                   ...
                    - i*(exchange/diffusion)*(SpinVector.Pauli*g - g*conj(SpinVector.Pauli))    ...
                    + 2i*(spinorbit.z + g*conj(spinorbit.z)*gt)*N*dg                            ...
                    + 2i*dg*Nt*(conj(spinorbit.z) + gt*spinorbit.z*g)                           ...
                    + 2*(spinorbit*g + g*conj(spinorbit))*Nt*(conj(spinorbit) + gt*spinorbit*g) ...
                    + (spinorbit^2*g - g*conj(spinorbit)^2);
            
            d2gt =  - 2*dgt*N*g*dgt                                                             ...
                    - 2i*(energy/diffusion)*gt                                                  ...
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
        
        
        function boundary()
        end
    end
end