% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-14
% Updated 2015-02-15

function A = RashbaDresselhaus(strength, angle)
    % This function returns an SU(2) vector field that describes the
    % a Rashba--Dresselhaus coupling in the xy-plane. The coupling
    % constants are given in polar coordinates, so that the Rashba constant
    % is strength*sin(angle), and the Dresselhaus one strength*cos(angle).
    %
    % Input:
    %   strength        Strength of the spin-orbit coupling;
    %                   increases both Rashba and Dresselhaus couplings.
    %   angle           Angle between coupled spin and momentum components;
    %                   rotates between Dresselhaus and Rashba couplings.
    % Output:
    %   A               3x2x2 SU(2) valued vector field that describes the
    %                   Rashba--Dresselhaus spin-orbit coupling above.

    % Define the Pauli vector and matrices
    Pauli = SpinVector([0,1; 1,0], [0,-1i; 1i,0], [1,0; 0,-1]);

    % Define the Rashba--Dresselhaus SU(2) field
    A = SpinVector( strength*(cos(angle)*Pauli.x + sin(angle)*Pauli.y), ...
                   -strength*(cos(angle)*Pauli.y + sin(angle)*Pauli.x), ...
                    0 );
end