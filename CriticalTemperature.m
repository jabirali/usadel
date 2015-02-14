function gap = CriticalTemperature(singlet, temperature, cutoff)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here\

    % BCS test:
    gap     = singlet
    singlet = @(energy) gap./sqrt(gap.^2 + energy.^2)

    kernel = @(energy) singlet(energy).*tanh(energy/(2.*temperature))
    gap    = integral(kernel, 0, cutoff)

end

