% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-15
% Updated 2015-02-15

function gap = SuperconductingGap(singlet, temperature, cutoff, scaling)
    % This function takes a the singlet component of the Green's function
    % in a superconductor and the temperature of the system, and calculates
    % the superconducting gap of the material.
    %
    % Input:
    %           scaling     The superconducting gap is proportional to this
    %                       scaling constant (which is usually called N₀λ)

    %singlet = @(energy) gap./sqrt(gap.^2 + energy.^2) x-dependent!

    int = @(x,energy) real(singlet(x,energy)) .* tanh(energy/(2.*temperature));
    gap = @(x) integral(int(x,energy), 0, cutoff);
end