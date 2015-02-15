% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-15
% Updated 2015-02-15


% Should I create a class containing JacobianS, JacobianF, etc?

function states = N(positions, states, energy, gap, exchange, spinorbit, temperature)
% Solve Usadel eq self-consistently
    addpath('BVP/');

    %options  = bvpset('AbsTol',1e-06,'RelTol',1e-06,'Nmax',2000);
    %solinit  = bvpinit(positions, states);
    %system   = bvp6c(@mat4ode,@mat4bc,solinit,options);
    %deeval
    
    for i=1:1000
        r = jacobianF(0,State(spinorbit.x,spinorbit.y,spinorbit.y,spinorbit.z).vectorize);
    end
    
    for n=1:length(states)
        solution(:,n) = states(n).vectorize;
    end
    
    for n=1:length(states)
        states(n) = State(solution(:,n));
    end
    % 
    %
    %
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subfunctions describing the behaviour of a superconductor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dydx = jacobianS(x,y)
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

    function boundaryS()
        % Calculates Kuprianov--Lukichev
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subfunctions describing the behaviour of a ferromagnet
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dydx = jacobianF(x,y)
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
                - 2i*energy*g                                                               ...
                - i*exchange*(SpinVector.Pauli*g - g*conj(SpinVector.Pauli))                ...
                + 2i*(spinorbit.z + g*conj(spinorbit.z)*gt)*N*dg                            ...
                + 2i*dg*Nt*(conj(spinorbit.z) + gt*spinorbit.z*g)                           ...
                + 2*(spinorbit*g + g*conj(spinorbit))*Nt*(conj(spinorbit) + gt*spinorbit*g) ...
                + (spinorbit^2*g - g*conj(spinorbit)^2);
        
        d2gt =  - 2*dgt*N*g*dgt                                                             ...
                - 2i*energy*gt                                                              ...
                + i*exchange*(conj(SpinVector.Pauli)*gt - gt*SpinVector.Pauli)              ...
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
end