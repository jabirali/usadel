function [singlet,tripletS,tripletL] = decomposition(green, exchange)
    % This function takes a Green's function matrix in spin space as input,
    % and splits the Green's function into a singlet component, a
    % short-range triplet component (along the exchange field), and a
    % long-range triplet component (perpendicular to the exchange field).
    % 
    % Input:    green       2x2 Green's function matrix
    %           exchange    3x1 exchange field vector
    %
    % Output:   singlet     1x1 singlet component
    %           tripletS    3x1 short-range triplet component
    %           tripletL    3x1 long-range triplet component
    
    % Normalize the exchange field
    unitvec = exchange/norm(exchange);

    % Singlet component (proportional to iσ^y)
    singlet  = (green(1,2) - green(2,1))/2;

    % Calculate the triplet component (proportional to [σ^x,σ^y,σ^z]iσ^y)
    triplet  = [(green(2,2) - green(1,1))/2;
                (green(1,1) + green(2,2))/2i;
                (green(1,2) + green(2,1))/2];
    
    % Project the triplet component along the exchange field to obtain the
    % short-range component, and the rest is then the long-range component
    tripletS = dot(triplet,unitvec) .* unitvec;
    tripletL = triplet - tripletS; 
end