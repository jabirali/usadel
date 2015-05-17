% This defines a data structure that describes the physical state of a
% metal for a given range of positions and energies. The purpose of this 
% class is mainly to be used as a base class for more interesting material
% classes, such as those that describe superconductors and ferromagnets,
% or to be used in conjunction with such materials in hybrid structures.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-23
% Updated 2015-05-07

classdef Metal < handle
    properties (GetAccess=public, SetAccess=public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define the internal variables of the Metal class
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Properties used to store the internal state of the system
        positions       = [];                % Vector of positions in the material
        energies        = [];                % Vector of energies for the material
        states          = State.empty(0,0);  % Riccati parameters and derivatives for each position and energy
        boundary_left   = State.empty(0);    % Riccati parameters and derivatives at the left boundary
        boundary_right  = State.empty(0);    % Riccati parameters and derivatives at the right boundary
        
        % Properties that determine the physical characteristics of the system
        interface_left  = inf;               % Interface parameter zeta at the left boundary
        interface_right = inf;               % Interface parameter zeta at the right boundary
        transparent     = false;             % Whether we should use transparent boundary conditions
        spinactive      = false;             % Whether we should use spin-active boundary conditions
        thouless        = 1;                 % Thouless energy of the system
        
        % Properties that describe spin-active interfaces
        magnetization_left  = [0,0,0];      % Magnetization vector at left interface  (unit vector)
        magnetization_right = [0,0,0];      % Magnetization vector at right interface (unit vector)
        polarization_left   = 0;            % Polarization at left interface  [-1,+1]
        polarization_right  = 0;            % Polarization at right interface [-1,+1]
        phaseshift_left     = 0;            % Spin-dependent phase-shift factor at left interface  [-inf,inf]
        phaseshift_right    = 0;            % Spin-dependent phase-shift factor at right interface [-inf,inf]
        
        % Internal properties that will be set and used during simulations
        coeff1  = {};                        % Coefficients in the differential equations for gamma
        coeff2  = {};                        % Coefficients in the differential equations for gamma~
        jc      = {};                        % Partially evaluated Jacobian functions
        bc      = {};                        % Partially evaluated boundary conditions
        
        % Properties that determine the simulation behavior
        error_abs = 1e-6;                    % Maximum absolute error when simulating
        error_rel = 1e-6;                    % Maximum relative error when simulating
        grid_size = 32768;                   % Maximum grid size to use in simulations
        
        % Properties that determine the behaviour of the program
        debug = true;                        % Whether to show intermediate results or not
        plot  = true;                        % Whether to plot intermediate results or not
        delay = 0;                           % How long to wait between program iterations
    end
    
    
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods that instantiate the object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = Metal(positions, energies, thouless)
            % This method constructs a Metal instance from a vector of
            % positions, a vector of energies, and the Thouless energy.
            
            % Set the internal properties to the provided values
            self.positions = positions;
            self.energies  = energies;
            self.thouless  = thouless;

            % Initialize the internal state of the metal
            self.states(length(self.positions), length(self.energies)) = State;
            self.init_superconductor;

            % Set the boundary conditions to vacuum states
            self.boundary_left(length(energies))  = State;
            self.boundary_right(length(energies)) = State;
        end

        function init_metal(self)
            % Initialize the internal state to a normal metal. This is useful
            % as an initial guess when simulating *weak* proximity effects.
            for i=1:length(self.positions)
               for j=1:length(self.energies)
                   self.states(i,j) = State;
               end
            end
        end
        
        function init_superconductor(self)
            % Initialize the internal state to a bulk superconductor with
            % superconducting gap 1. This is useful as an initial guess
            % when simulating *strong* proximity effects in energy units 
            % where the zero-temperature gap is normalized to unity.
            for i=1:length(self.positions)
               for j=1:length(self.energies)
                   self.states(i,j) = Superconductor.Bulk(self.energies(j), 1, 0);
               end
            end
        end
            


        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods that update the internal state of the object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function update_coeff(self)
            % This function updates the functions and coefficients passed to
            % the functions 'jacobian' and 'boundary' when solving equations.
            
            % Partially evaluate the Jacobian matrix and boundary conditions
            % for the different material energies, and store the resulting 
            % anonymous functions in an array. These functions are passed on
            % to bvp6c when solving the equations.
            self.jc = {};
            self.bc = {};
            for m=1:length(self.energies)
                self.jc{m} = @(x,y) self.jacobian(self,x,y,self.energies(m));
                if self.transparent
                    % If 'transparent' is true, use transparent b.c.
                    self.bc{m} = @(a,b) self.boundary_transparent(self,a,b,self.energies(m));
                elseif self.spinactive
                    % If 'spinactive' is true, use spin-active b.c.
                    self.bc{m} = @(a,b) self.boundary_spinactive(self,a,b,self.energies(m));
                else
                    % Else, use standard Kuprianov-Lukichev b.c. instead
                    self.bc{m} = @(a,b) self.boundary(self,a,b,self.energies(m));
                end
            end
        end
        
        function update_boundary_left(self, other)
            % This function updates the boundary condition to the left
            % based on the current state of an adjoining material.
            self.boundary_left(:) = other.states(end,:);
        end

        function update_boundary_right(self, other)
            % This function updates the boundary condition to the right
            % based on the current state of an adjoining material.
            self.boundary_right(:) = other.states(1,:);
        end

        function update_state(self)
            % This function solves the Usadel equation numerically for the
            % given position and energy range, and updates the current
            % estimate for the state of the system.
            
            % Set the accuracy of the numerical solution
            options = bvpset('AbsTol',self.error_abs,'RelTol',self.error_rel,'Nmax',self.grid_size);

            for m=1:length(self.energies)
                % Progress information
                self.print('[ %3.f / %3.f ]  E = %2.4f ', m, length(self.energies), self.energies(m));
                
                % Vectorize the current state of the system for the given
                % energy, and use it as an initial guess for the solution
                current = zeros(16,length(self.positions));
                for n=1:length(self.positions)
                    current(:,n) = self.states(n,m).vectorize;
                end
                initial = bvpinit(self.positions', current);
                
                % Solve the differential equation, and evaluate the
                % solution on the position vector of the metal
                solution = deval(bvp6c(self.jc{m},self.bc{m},initial,options), self.positions);
                
                % Update the current state of the system based on the solution
                for n=1:length(self.positions)
                    self.states(n,m) = State(solution(:,n));
                end
                
                % Time delay between iterations (reduces load on the system)
                pause(self.delay);
            end
        end
        
        
        function update(self)
            % This function updates the internal state of the
            % metal object by calling the other update methods.
            
            % Update the state
            self.update_coeff;
            self.update_state;
                
            % Plot the current density of states (if 'plot' is set to true)
            if self.plot
                self.plot_dos_surf;
            end
        end
    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define miscellaneous helper methods for the class
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function print(self,varargin)
            % This function is used to print messages, such as debug info
            % and progress, if the 'debug' flag is set to 'true'. The
            % message is preceded by the name of the class, which lets you
            % distinguish the output of 'Metal' instances from the output
            % of instances of daughter classes.
            
            if self.debug
                fprintf(':: %s: %s\n', class(self), sprintf(varargin{:}));
            end
        end

        function index = position_index(self, position)
            % Returns the vector index corresponding to a given energy value
            index = find(abs(self.positions-position) < 1e-8, 1, 'first');
        end
        
        function index = energy_index(self, energy)
            % Returns the vector index corresponding to a given energy value
            index = find(abs(self.energies-energy) < 1e-8, 1, 'first');
        end

        function backup = backup_save(self)
            % Returns a backup of the state of the metal
            backup = self.states(:,:);
        end
        
        function backup_load(self, backup)
            % Restores the state of the metal from a backup
            % NB: Remember to run 'update' after a call to 'backup_load' to
            %     make sure that the rest of the properties are consistent!
            self.states(:,:) = backup;
        end
        
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define methods for printing and plotting the internal state
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function plot_dos_left(self)
            % Calculate the *left* density of states for the system.
            % Note: we assume that we only have data for positive energies,
            %       and that the negative energy region is symmetric.

            N = length(self.positions);
            M = length(self.energies);
            dos = zeros(2*M-1);
            erg = zeros(2*M-1);
           
            for m=1:M
                erg(M-m+1) = -self.energies(m);
                dos(M-m+1) = self.states(1,m).eval_ldos;
                
                erg(M+m-1) = self.energies(m);
                dos(M+m-1) = dos(M-m+1);
            end
            
            % Plot the results
            plot(erg, dos);
            xlabel('\epsilon/\Delta_0')
            ylabel('D(\epsilon)')
            set(gca, 'XTick', -20:20);
            set(gca, 'YTick', -20:20);
            set(gca, 'XLim',  [-2 2]);
            set(gca, 'YLim',  [ 0 2]);
        end
        
        function plot_dos_right(self)
            % Calculate the *right* density of states for the system.
            % Note: we assume that we only have data for positive energies,
            %       and that the negative energy region is symmetric.

            N = length(self.positions);
            M = length(self.energies);
            dos = zeros(2*M-1);
            erg = zeros(2*M-1);
            for m=1:M
                erg(M-m+1) = -self.energies(m);
                dos(M-m+1) = self.states(end,m).eval_ldos;
                
                erg(M+m-1) = self.energies(m);
                dos(M+m-1) = dos(M-m+1);
            end
            
            % Plot the results
            plot(erg, dos);
            xlabel('\epsilon/\Delta_0')
            ylabel('D(\epsilon)')
            set(gca, 'XTick', -20:20);
            set(gca, 'YTick', -20:20);
            set(gca, 'XLim',  [-2 2]);
            set(gca, 'YLim',  [0 2]);
        end
        
        function plot_dos_center(self)
            % Calculate the *central* density of states for the system.
            % Note: we assume that we only have data for positive energies,
            %       and that the negative energy region is symmetric.

            N = length(self.positions);
            M = length(self.energies);
            dos = zeros(2*M-1);
            erg = zeros(2*M-1);
            pos = floor(length(self.positions)/2);
            for m=1:M
                erg(M-m+1) = -self.energies(m);
                dos(M-m+1) = self.states(pos,m).eval_ldos;
                
                erg(M+m-1) = self.energies(m);
                dos(M+m-1) = dos(M-m+1);
            end
            
            % Plot the results
            plot(erg, dos);
            xlabel('\epsilon/\Delta_0')
            ylabel('D(\epsilon)')
            set(gca, 'XTick', -20:20);
            set(gca, 'YTick', -20:20);
            set(gca, 'XLim',  [-2 2]);
            set(gca, 'YLim',  [0 2]);
        end
        
        function plot_dos_surf(self)
            % Calculate the density of states throughout the system.
            % Note: we assume that we only have data for positive energies,
            %       and that the negative energy region is symmetric.
            
            N = length(self.positions);
            M = length(self.energies);
            dos = zeros(N, 2*M-1);
            for n=1:N
                dos(n,M) = self.states(n,1).eval_ldos;
                for m=2:M
                    dos(n,M-m+1) = self.states(n,m).eval_ldos;
                    dos(n,M+m-1) = dos(n,M-m+1);
                end
            end
                        
            % Plot the results
            surf([fliplr(-self.energies) self.energies(2:end)], self.positions, dos, 'EdgeColor', 'none');
            shading('interp');
            colormap(parula(256));
            caxis([0 2]);
            view(7.5,30);
            
            set(gca, 'XTick', -20:20);
            set(gca, 'YTick', -20:20);
            set(gca, 'ZTick', -20:20);
            set(gca, 'XLim',  [-2 2 ]);
            set(gca, 'ZLim',  [ 0 2 ]);
        end
        
        function plot_dist(self)
            % Calculate the singlet and triplet distributions.
            %
            % NB: This implementation only *adds* the contribution from
            %     every energy, and does not perform a proper integral!

            singlet = zeros(length(self.positions), 1);
            triplet = zeros(length(self.positions), 1);
            for m=1:length(self.energies)
                for n=1:length(self.positions)
                    singlet(n) = singlet(n) + norm(self.states(n,m).singlet);
                    triplet(n) = triplet(n) + norm(self.states(n,m).triplet);
                end
            end
            
            % Plot cubic interpolations of the results
            positions = linspace(self.positions(1), self.positions(end), 100);
            plot(positions, pchip(self.positions, singlet, positions), ...
                 positions, pchip(self.positions, triplet, positions));
            xlabel('Relative position');
            ylabel('Distribution');
            legend('Singlet', 'Triplet');
        end
    end
    
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define static methods (available without object instantiation)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function dydx = jacobian(self, x, y, energy)
            % This function takes a Metal object 'self', the position 'x',
            % the current state vector 'y', and an energy as inputs, and
            % calculates the Jacobian of the system. This is performed
            % using the Riccati parametrized Usadel equations.
            
            % Extract the Riccati parameters and their derivatives
            [g,dg,gt,dgt] = State.unpack(y);
            
            % Calculate the normalization matrices
            N  = inv( eye(2) - g*gt );
            Nt = inv( eye(2) - gt*g );
            
            % Calculate the second derivatives of the Riccati parameters
            % according to the Usadel equation in the metal
            d2g  = -2 * dg*Nt*gt*dg - (2i/self.thouless) * (energy+1e-3i)*g;
            d2gt = -2 * dgt*N*g*dgt - (2i/self.thouless) * (energy+1e-3i)*gt;
            
            % Pack the results into a state vector
            dydx = State.pack(dg,d2g,dgt,d2gt);
        end
        
        function residue = boundary(self, y1, y2, energy)
            % This function takes a Metal object 'self', the position 'x',
            % the current state vector 'y', and an energy as inputs, and
            % calculates the Kuprianov-Lukichev boundary conditions. This
            % function will be used as b.c. when 'transparent' is false.
            
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
            dg1  = dg1  - (1/self.interface_left)*( eye(2) - g1*gt0 )*N0*(  g1  - g0  );
            dgt1 = dgt1 - (1/self.interface_left)*( eye(2) - gt1*g0 )*Nt0*( gt1 - gt0 );
            
            dg2  = dg2  - (1/self.interface_right)*( eye(2) - g2*gt3 )*N3*(  g3  - g2  );
            dgt2 = dgt2 - (1/self.interface_right)*( eye(2) - gt2*g3 )*Nt3*( gt3 - gt2 );
            
            % Vectorize the results of the calculations, and return it
            residue = State.pack(dg1,dgt1,dg2,dgt2);
        end

        function residue = boundary_spinactive(self, y1, y2, energy)
            % This function takes a Metal object 'self', the position 'x',
            % the current state vector 'y', and an energy as inputs, and
            % calculates the spin-active version of the Kuprianov-Lukichev
            % boundary conditions. This function will be used as b.c. when 
            % the property 'spinactive' is set to 'true'.
            
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
            
            % Calculate the Green's functions
            I = eye(2);
            
            G0  = (I - g0*gt0) \ (I + g0*gt0);
            Gt0 = (I - gt0*g0) \ (I + gt0*g0);
            F0  = (I - g0*gt0) \ (2*g0);
            Ft0 = (I - gt0*g0) \ (2*gt0);
            
            G1  = (I - g1*gt1) \ (I + g1*gt1);
            Gt1 = (I - gt1*g1) \ (I + gt1*g1);
            F1  = (I - g1*gt1) \ (2*g1);
            Ft1 = (I - gt1*g1) \ (2*gt1);

            G2  = (I - g2*gt2) \ (I + g2*gt2);
            Gt2 = (I - gt2*g2) \ (I + gt2*g2);
            F2  = (I - g2*gt2) \ (2*g2);
            Ft2 = (I - gt2*g2) \ (2*gt2);

            G3  = (I - g3*gt3) \ (I + g3*gt3);
            Gt3 = (I - gt3*g3) \ (I + gt3*g3);
            F3  = (I - g3*gt3) \ (2*g3);
            Ft3 = (I - gt3*g3) \ (2*gt3);
            
            % Calculate the interface parameters
            LS = self.magnetization_left * SpinVector.Pauli;                                  % L: m·σ
            LT = self.magnetization_left * conj(SpinVector.Pauli);                            % L: m·σ*
            LM = self.polarization_left/(1 + sqrt(1-self.polarization_left^2));               % L: µ
            LK = (1-sqrt(1-self.polarization_left^2))/(1+sqrt(1-self.polarization_left^2));   % L: κ
            LL = 1i*self.phaseshift_left;                                                     % L: iλ
            
            RS = self.magnetization_right * SpinVector.Pauli;                                 % R: m·σ
            RT = self.magnetization_right * conj(SpinVector.Pauli);                           % R: m·σ*
            RM = self.polarization_right/(1 + sqrt(1-self.polarization_right^2));             % R: µ
            RK = (1-sqrt(1-self.polarization_right^2))/(1+sqrt(1-self.polarization_right^2)); % R: κ
            RL = 1i*self.phaseshift_right;                                                    % R: iλ

            % Calculate the left interface matrices
            L1  = (G1*G0 - F1*Ft0)*(I+LM*LS) - (I+LM*LS)*(G0*G1 - F0*Ft1)                        ...
                + (G1*LS*G0 - F1*LT*Ft0)*(LM + LK*LS) - (LM + LK*LS)*(G0*LS*G1 - F0*LT*Ft1)      ...
                + LL*(G1*LS - LS*G1);
           
            Lt1 = (Gt1*Gt0 - Ft1*F0)*(I+LM*LT) - (I+LM*LT)*(Gt0*Gt1 - Ft0*F1)                    ...
                + (Gt1*LT*Gt0 - Ft1*LS*F0)*(LM + LK*LT) - (LM + LK*LT)*(Gt0*LT*Gt1 - Ft0*LS*F1)  ...
                - LL*(Gt1*LT - LT*Gt1);
            
            L2  = (G1*F0 - F1*Gt0)*(I+LM*LT) - (I+LM*LS)*(G0*F1 - F0*Gt1)                        ...
                + (G1*LS*F0 - F1*LT*Gt0)*(LM+LK*LT) - (LM+LK*LS)*(G0*LS*F1 - F0*LT*Gt1)          ...
                + LL*(F1*LT - LS*F1);

            Lt2 = (Gt1*Ft0 - Ft1*G0)*(I+LM*LS) - (I+LM*LT)*(Gt0*Ft1 - Ft0*G1)                    ...
                + (Gt1*LT*Ft0 - Ft1*LS*G0)*(LM+LK*LS) - (LM+LK*LT)*(Gt0*LT*Ft1 - Ft0*LS*G1)      ...
                - LL*(Ft1*LS - LT*Ft1);            
            
            % Calculate the right interface matrices
            R1  = (G2*G3 - F2*Ft3)*(I+RM*RS) - (I+RM*RS)*(G3*G2 - F3*Ft2)                        ...
                + (G2*RS*G3 - F2*RT*Ft3)*(RM + RK*RS) - (RM + RK*RS)*(G3*RS*G2 - F3*RT*Ft2)      ...
                + RL*(G2*RS - RS*G2);
           
            Rt1 = (Gt2*Gt3 - Ft2*F3)*(I+RM*RT) - (I+RM*RT)*(Gt3*Gt2 - Ft3*F2)                    ...
                + (Gt2*RT*Gt3 - Ft2*RS*F3)*(RM + RK*RT) - (RM + RK*RT)*(Gt3*RT*Gt2 - Ft3*RS*F2)  ...
                - RL*(Gt2*RT - RT*Gt2);
            
            R2  = (G2*F3 - F2*Gt3)*(I+RM*RT) - (I+RM*RS)*(G3*F2 - F3*Gt2)                        ...
                + (G2*RS*F3 - F2*RT*Gt3)*(RM+RK*RT) - (RM+RK*RS)*(G3*RS*F2 - F3*RT*Gt2)          ...
                + RL*(F2*RT - RS*F2);

            Rt2 = (Gt2*Ft3 - Ft2*G3)*(I+RM*RS) - (I+RM*RT)*(Gt3*Ft2 - Ft3*G2)                    ...
                + (Gt2*RT*Ft3 - Ft2*RS*G3)*(RM+RK*RS) - (RM+RK*RT)*(Gt3*RT*Ft2 - Ft3*RS*G2)      ...
                - RL*(Ft2*RS - RT*Ft2);
            
            % Calculate the deviation from the Kuprianov--Lukichev b.c.
            % with spin-active interface terms from Machon et al.
            dg1  = dg1  + (0.25/self.interface_left)*(I - g1*gt1)*(L2  - L1*g1);
            dgt1 = dgt1 + (0.25/self.interface_left)*(I - gt1*g1)*(Lt2 - Lt1*gt1);
            
            dg2  = dg2  - (0.25/self.interface_right)*(I - g2*gt2)*(R2 - R1*g2);
            dgt2 = dgt2 - (0.25/self.interface_right)*(I - gt2*g2)*(Rt2 - Rt1*gt2);
            
            % Vectorize the results of the calculations, and return it
            residue = State.pack(dg1,dgt1,dg2,dgt2);
        end
        
        function residue = boundary_transparent(self, y1, y2, energy)
            % This function takes a Metal object 'self', the position 'x',
            % the current state vector 'y', and an energy as inputs, and
            % calculates the transparent boundary conditions. This function
            % will be used as b.c. when 'transparent' is set to true. Note
            % that if the interface parameters are infinite, we assume a
            % vacuum interface, and therefore use the derivative boundary
            % condition dg = 0 instead of the transparency condition.
            
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
            
            % Calculate the deviation from the boundary conditions. We use
            % the derivative boundary condition dg = 0 when there is an
            % infinite interface parameter (i.e. vacuum interface), and the
            % transparent boundary condition g1 = g2/interface_param when
            % there is a finite interface parameter.

            if isinf(self.interface_left)
                r1  = dg1;
                rt1 = dgt1;
            else
                r1  = g1  - g0/self.interface_left;
                rt1 = gt1 - gt0/self.interface_left;
            end
            
            if isinf(self.interface_right)
                r2  = dg2;
                rt2 = dgt2;
            else
                r2  = g2  - g3/self.interface_right;
                rt2 = gt2 - gt3/self.interface_right;
            end
            
            % Vectorize the results of the calculations, and return it
            residue = State.pack(r1,rt1,r2,rt2);
        end
    end
end