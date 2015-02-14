% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Based on a similar program by Dr. Sol Jacobsen
% Created 2015-02-15
% Updated 2015-02-15

function [gg1,gg2,gt,dgt] = UsadelF(state, exchange, spinorbit)
    % This function
    
    % Make sure the BVP package files are in the current Matlab search path
    addpath('BVP/')
    
    % Initialize the boundary value problem solve with an initial guess
    % given by the function 'mat4init', a from
    % 'mat4init' 'mat4initThe first command sets
    % the error tolerances, the second defines the initial guess for the
    % solution using the function 'mat4init', and the final one solves the
    % differential equation using the Jacobian 'mat4ode' and boundary cond


    
    options  = bvpset('AbsTol',1e-06,'RelTol',1e-06,'Nmax',2000);            % Error tolerances
    solinit  = bvpinit(linspace(0,1,100),@mat4init);                         % Initial guess given by 'mat4init'
    system   = bvp6c(@mat4ode,@mat4bc,state,options);                          % Solve the differential equation using the Jacobian 'mat4ode',
                                                                            % boundary conditions 'mat4bc', and initial guess 'solinit'


    xs = linspace(0,1);%We want 100 solutions (sol) for different x between 0 and 1.
    ys = deval(system,xs);%deval returns solutions sol to a differential equation, using the bvp6c package above
        
    % Extract the gamma matrices and their derivatives drom the state vector
    gs   = reshape(ys( 1: 4), 2, 2)';
    dgs  = reshape(ys( 5: 8), 2, 2)';
    gts  = reshape(ys( 9:12), 2, 2)';
    dgts = reshape(ys(13:16), 2, 2)';
end


function dydx = mat4ode(x,y)
    % Extract the gamma matrices and their derivatives drom the state vector
    g   = reshape(y( 1: 4), 2, 2)';
    dg  = reshape(y( 5: 8), 2, 2)';
    gt  = reshape(y( 9:12), 2, 2)';
    dgt = reshape(y(13:16), 2, 2)';

    % Calculate the normalization matrices
    N  = inv( eye(2) - g*gt );
    Nt = inv( eye(2) - gt*g );


% We will now write the Usadel equations for gamma, tilde{gamma} in the
% form partial_x^2 \gamma = tot. The tot-matrix then includes first
% order-derivatives and the terms depending on E, h, and Delta. 

M = dg*Nt*gt*dg;
Mt = dgt*N*g*dgt;

spin  = 2i*E*g/Ethf  + i*h*Pauli*g/Ethf  - i*g*h*conj(Pauli)/Ethf;
spint = 2i*E*gt/Ethf + i*gt*h*Pauli/Ethf - i*h*conj(Pauli)*gt/Ethf;

% Here define the additional terms due to the INTRINSIC spin-orbit
% coupling, again divide by Thouless energy in ferromagnet for
% normalisation...

%Define some component cell arrays for the second SO term:
ant1={A{1}*g,A{2}*g,A{3}*g};
ant2={g*conj(A{1}),g*conj(A{2}),g*conj(A{3})};
ant={ant1{1}+ant2{1},ant1{2}+ant2{2},ant1{3}+ant2{3}};
part1={conj(A{1})+gt*A{1}*g, conj(A{2})+gt*A{2}*g, conj(A{3})+gt*A{3}*g};
part2={Nt*part1{1},Nt*part1{2},Nt*part1{3}};

ant1t={conj(A{1})*gt,conj(A{2})*gt,conj(A{3})*gt};
ant2t={gt*A{1},gt*A{2},gt*A{3}};
antt={ant1t{1}+ant2t{1},ant1t{2}+ant2t{2},ant1t{3}+ant2t{3}};
part1t={A{1}+g*conj(A{1})*gt, A{2}+g*conj(A{2})*gt, A{3}+g*conj(A{3})*gt};
part2t={N*part1t{1},N*part1t{2},N*part1t{3}};

%Remember I have the "wrong" normalisation here as I divide through by Ethf
orb = A^2 * g - g * conj(A)^2 + 2*(ant{1}*part2{1}+ant{2}*part2{2}+ant{3}*part2{3}) + 2i*(dg*Nt*(conj(A{3})+gt*A{3}*g)+(A{3}+g*conj(A{3})*gt)*N*dg))/Ethf; 
orbt = conj(A)^2 * gt - gt * A^2 + 2*(antt{1}*part2t{1}+antt{2}*part2t{2}+antt{3}*part2t{3}) - 2i*(dgt*N*(A{3}+g*conj(A{3})*gt)+(conj(A{3})+gt*A{3}*g)*Nt*dgt))/Ethf;

tot = -2*M - spin + orb;
tott = -2*Mt-spint + orbt;

tot11 = tot(1,1);
tot12 = tot(1,2);
tot21 = tot(2,1);
tot22 = tot(2,2);
tott11 = tott(1,1);
tott12 = tott(1,2);
tott21 = tott(2,1);
tott22 = tott(2,2);

dydx = [dg rhs ]

dydx = [dg1; tot11;
    dg2; tot12;
    dg3; tot21;
    dg4; tot22;
    dgt1; tott11;
    dgt2; tott12;
    dgt3; tott21;
    dgt4; tott22];

end

% =======================================================================
function res = mat4bc(ya,yb)

% Vil bruke verdien av Green's function i S  ved interface i
% grensebetingelsen
gs1 = ggs1(jj);
gs2 = ggs2(jj); 
gs3 = ggs3(jj);
gs4 = ggs4(jj); 

gts1 = ggts1(jj);
gts2 = ggts2(jj); 
gts3 = ggts3(jj);
gts4 = ggts4(jj); 

% At x=0
g1_x0 = ya(1);
g2_x0 = ya(3);
g3_x0 = ya(5);
g4_x0 = ya(7);
gt1_x0 = ya(9);
gt2_x0 = ya(11);
gt3_x0 = ya(13);
gt4_x0 = ya(15);

dg1_x0 = ya(2);
dg2_x0 = ya(4);
dg3_x0 = ya(6);
dg4_x0 = ya(8);
dgt1_x0 = ya(10);
dgt2_x0 = ya(12);
dgt3_x0 = ya(14);
dgt4_x0 = ya(16);

g_x0 = [g1_x0,g2_x0; g3_x0,g4_x0];
gt_x0 = [gt1_x0,gt2_x0; gt3_x0,gt4_x0];

gs = [gs1,gs2;gs3,gs4];
gts = [gts1,gts2;gts3,gts4];


% At x=1
g1_xd = yb(1);
g2_xd = yb(3);
g3_xd = yb(5);
g4_xd = yb(7);
gt1_xd = yb(9);
gt2_xd = yb(11);
gt3_xd = yb(13);
gt4_xd = yb(15);

dg1_xd = yb(2);
dg2_xd = yb(4);
dg3_xd = yb(6);
dg4_xd = yb(8);
dgt1_xd = yb(10);
dgt2_xd = yb(12);
dgt3_xd = yb(14);
dgt4_xd = yb(16);

g_xd = [g1_xd,g2_xd; g3_xd,g4_xd];
gt_xd = [gt1_xd,gt2_xd; gt3_xd,gt4_xd];

id = [1,0;0,1];

Nj = inv(id-g_x0*gt_x0);
Ntj = inv(id-gt_x0*g_x0);

Ns = inv(id-gs*gts);
Nts = inv(id-gts*gs);

%Note that in the BCs I take explicitly the Az component coded here for
%future-proofing the code if we want to test noncentrosym. materials but
%this component is acutally zero in my case so you can just take it out if
%you want.

% Tot is now the expression for my derivative of gamma_2, ie the
% ferromagnetic gamma, at x=0. Tott is the tilded version. 
Tot=((id-g_x0*gts)*Ns*(g_x0-gs))/(zeta)+i*A{3}*g_x0+i*g_x0*conj(A{3});
Tott=((id-gt_x0*gs)*Nts*(gt_x0-gts))/(zeta)-i*conj(A{3})*gt_x0-i*gt_x0*A{3};

Tot2= i*A{3}*g_xd+i*g_xd*conj(A{3});
Tott2=-i*conj(A{3})*gt_xd-i*gt_xd*A{3};

% Boundary conditions summarized. Using KL b.c. at the S/F interface
% ie we have derivative of first element of ferromagnetic gamma at x=0 is given by the
% first element of Tot. At x=L(=1) the derivative zero since we have a
% vacuum there
res = [dg1_x0 - Tot(1,1);
    dg2_x0 - Tot(1,2);
    dg3_x0 - Tot(2,1);
    dg4_x0 - Tot(2,2);
    dgt1_x0 - Tott(1,1);
    dgt2_x0 - Tott(1,2);
    dgt3_x0 - Tott(2,1);
    dgt4_x0 - Tott(2,2);
    dg1_xd - Tot2(1,1);
    dg2_xd - Tot2(1,2);
    dg3_xd - Tot2(2,1);
    dg4_xd - Tot2(2,2);
    dgt1_xd - Tott2(1,1);
    dgt2_xd - Tott2(1,2);
    dgt3_xd - Tott2(2,1);
    dgt4_xd - Tott2(2,2)];z