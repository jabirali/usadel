% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-15
% Updated 2015-02-15

function gap = SuperconductingGap(states, energies, temperature, cutoff, scaling)
    % This function takes a the singlet component of the Green's function
    % in a superconductor and the temperature of the system, and calculates
    % the superconducting gap of the material.
    %
    % Input:    singlet     
    %           temperature The temperature of the material
    %           cutoff      The Debye frequency cutoff of the material
    %           scaling     The superconducting gap is proportional to this
    %                       scaling constant (which is usually called N₀λ)

    %singlet = @(energy) gap./sqrt(gap.^2 + energy.^2) x-dependent!

    
    % Calculate and collect the singlet components of the states
    singlets = zeros(1,numel(states));
    for n=1:length(singlets)
        singlets(n) = states(n).singlet;
    end

    % Create a cubic interpolation of the numerical data above, multiplied
    % by the tanh(ε/2T) kernel in the equation for the superconducting gap
    kernel = @(E) pchip(energies, singlets, E) .* tanh(E./(2*temperature));
    
    % Perform a numerical integration of the interpolation up to the cutoff
    gap    = scaling * integral(kernel, 0, cutoff);
    
    %int = @(x,energy) real(singlet(x,energy)) .* tanh(energy/(2.*temperature));
    %gap = @(x) integral(int(x,energy), 0, cutoff);
end