% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-14
% Updated 2015-02-15

function gap = SuperconductingGap(states, energies, temperature, cutoff, scaling)
    % This function takes a the singlet component of the Green's function
    % in a superconductor and the temperature of the system, and calculates
    % the superconducting gap of the material.
    %
    % Input:   
    %   states      This should be a vector containing either State objects
    %               (numerical solution of the Usadel eq) or complex numbers
    %               (the singlet components of the solutions).
    %   energies    This should be a vector of energies which corresponds
    %               to the states above.
    %   temperature The temperature of the material.
    %   cutoff      The Debye frequency cutoff of the material.
    %   scaling     The superconducting gap is proportional to this
    %               scaling constant (which is usually written N₀λ).
    % 
    % Output:  
    %   gap         The superconducting gap (in general a complex number)
    
    % Extract the singlet components from the states, if necessary
    singlets = zeros(1,numel(states));
    if isa(states, 'State')
        for n=1:length(singlets)
            singlets(n) = states(n).singlet;
        end
    else
        singlets = states;
    end

    % Create a cubic interpolation of the numerical data above, multiplied
    % by the tanh(ε/2T) kernel in the equation for the superconducting gap
    kernel = @(E) pchip(energies, singlets, E) .* tanh(E./(2*temperature));
    
    % Perform a numerical integration of the interpolation up to the cutoff
    gap    = scaling * integral(kernel, 0, cutoff);
    
    %int = @(x,energy) real(singlet(x,energy)) .* tanh(energy/(2.*temperature));
    %gap = @(x) integral(int(x,energy), 0, cutoff);
end