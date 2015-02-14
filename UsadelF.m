function [gg1,gg2,gt,dgt] = UsadelF(exchange,spinorbit)
    % This function
    
    % Make sure the BVP package files are in the current Matlab search path
    addpath('BVP/')
    
    % Initialize the boundary value problem solver
    % The arguments to 'bvp6c' 
    options = bvpset('AbsTol',1e-06,'RelTol',1e-06,'Nmax',2000);
    solinit = bvpinit(linspace(0,1,100),@mat4init);
    sol = bvp6c(@mat4ode,@mat4bc,solinit,options);

%BVP is boundary value solver (installed separately)
%takes as its input functions mat4ode etx defined below. Function handle
%@ lets matlab know to expect a function as the argument here. The inputs
%are therefore the ode, boundary conditions and initial guess. Options sets
%some parameters for the calulational accuracy and number of steps.

xint = linspace(0,1);%We want 100 solutions (sol) for different x between 0 and 1.
Sxint = deval(sol,xint);%deval returns solutions sol to a differential equation, using the bvp6c package above

gg1 = Sxint(1,:);%ie here have solution 1 to diff eqn for all x values
gg2 = Sxint(3,:);%here have solution 3 for all x
gg3 = Sxint(5,:);
gg4 = Sxint(7,:);
ggt1 = Sxint(9,:);
ggt2 = Sxint(11,:);
ggt3 = Sxint(13,:);
ggt4 = Sxint(15,:);

dgg1 = Sxint(2,:);%the even solutions give the differential of gamma, which are also solutions of a bvp
dgg2 = Sxint(4,:);% the inputs in mat4ode and mat4bc thus include as evens the corresponding expressions for  
dgg3 = Sxint(6,:);%gamma, \tilde(gamma) and their derivatives
dgg4 = Sxint(8,:);
dggt1 = Sxint(10,:);
dggt2 = Sxint(12,:);
dggt3 = Sxint(14,:);
dggt4 = Sxint(16,:);

% ===================================================================
function dydx = mat4ode(x,y)

global Ethf
global E
global hx
global hy
global hz
global sigmax
global sigmay
global sigmaz
global A

% Define the elements of the gamma, tilde{gamma} matrices
% and their derivatives
g1 = y(1);
g2 = y(3);
g3 = y(5);
g4 = y(7);
gt1 = y(9);
gt2 = y(11);
gt3 = y(13);
gt4 = y(15);

dg1 = y(2);
dg2 = y(4);
dg3 = y(6);
dg4 = y(8);
dgt1 = y(10);
dgt2 = y(12);
dgt3 = y(14);
dgt4 = y(16);

% gamma, tilde{gamma}
g = [g1,g2;g3,g4];
gt = [gt1,gt2;gt3,gt4];

% Their derivatives
dg = [dg1, dg2; dg3, dg4];
dgt = [dgt1, dgt2; dgt3, dgt4];

% The elements of the normalization matrix N 
denom = (g1*gt1 + g2*gt3 + g3*gt2 + g4*gt4 - g1*g4*gt1*gt4 + g1*g4*gt2*gt3 + g2*g3*gt1*gt4 - g2*g3*gt2*gt3 - 1)^(-1);
N1 = denom*(g3*gt2 + g4*gt4 -1);
N2 = denom*(-g1*gt2 - g2*gt4);
N3 = denom*(-g3*gt1 - g4*gt3);
N4 = denom*(g1*gt1 + g2*gt3 - 1);

% The elements of the normalization matrix tilde{N}
denomt = (g1*gt1 + g2*gt3 + g3*gt2 + g4*gt4 - g1*g4*gt1*gt4 + g1*g4*gt2*gt3 + g2*g3*gt1*gt4 - g2*g3*gt2*gt3 - 1)^(-1);
Nt1 = denomt*(gt3*g2 + gt4*g4 - 1);
Nt2 = denomt*(-gt1*g2 - gt2*g4);
Nt3 = denomt*(-gt3*g1 - gt4*g3);
Nt4 = denomt*(gt1*g1 + gt2*g3 - 1);

N = [N1, N2; N3, N4];
Nt = [Nt1, Nt2; Nt3, Nt4];

% We will now write the Usadel equations for gamma, tilde{gamma} in the
% form partial_x^2 \gamma = tot. The tot-matrix then includes first
% order-derivatives and the terms depending on E, h, and Delta. 

M = dg*Nt*gt*dg;
Mt = dgt*N*g*dgt;

spin = 2*1i*E*g/Ethf + 1i*(hx*sigmax + hy*sigmay + hz*sigmaz)*g/Ethf - 1i*g*(hx*sigmax + hy*conj(sigmay) + hz*sigmaz)/Ethf;
spint = 2*1i*E*gt/Ethf + 1i*gt*(hx*sigmax + hy*sigmay + hz*sigmaz)/Ethf - 1i*(hx*sigmax + hy*conj(sigmay) + hz*sigmaz)*gt/Ethf;

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
orb = ((A{1}*A{1}+A{2}*A{2}+A{3}*A{3})*g-g*(conj(A{1})*conj(A{1})+conj(A{2})*conj(A{2})+conj(A{3})*conj(A{3})) + 2*(ant{1}*part2{1}+ant{2}*part2{2}+ant{3}*part2{3}) + 2i*(dg*Nt*(conj(A{3})+gt*A{3}*g)+(A{3}+g*conj(A{3})*gt)*N*dg))/Ethf; 
orbt = ((conj(A{1})*conj(A{1})+conj(A{2})*conj(A{2})+conj(A{3})*conj(A{3}))*gt-gt*(A{1}*A{1}+A{2}*A{2}+A{3}*A{3}) + 2*(antt{1}*part2t{1}+antt{2}*part2t{2}+antt{3}*part2t{3}) - 2i*(dgt*N*(A{3}+g*conj(A{3})*gt)+(conj(A{3})+gt*A{3}*g)*Nt*dgt))/Ethf;

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

dydx = [dg1; tot11;
    dg2; tot12;
    dg3; tot21;
    dg4; tot22;
    dgt1; tott11;
    dgt2; tott12;
    dgt3; tott21;
    dgt4; tott22];

% =======================================================================
function res = mat4bc(ya,yb)

global ggs1;
global ggs2;
global ggs3;
global ggs4;

global ggts1;
global ggts2;
global ggts3;
global ggts4;

global zeta
global jj

global A

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
    dgt4_xd - Tott2(2,2)];

% =========================================================================



function yinit = mat4init(x)

global jj

global ggf1_init;
global dggf1_init;
global ggf2_init;
global dggf2_init;
global ggf3_init;
global dggf3_init;
global ggf4_init;
global dggf4_init;

global ggtf1_init;
global dggtf1_init;
global ggtf2_init;
global dggtf2_init;
global ggtf3_init;
global dggtf3_init;
global ggtf4_init;
global dggtf4_init

k = round(100*x);

if k==0
    k=1;
end

% This is the guess for the y-vector, i.e.
% y=[gamma1,gamma1',gamma2,gamma2'...]

yinit = [ggf1_init(jj,k); 
         dggf1_init(jj,k); 
         ggf2_init(jj,k); 
         dggf2_init(jj,k); 
         ggf3_init(jj,k);  
         dggf3_init(jj,k); 
         ggf4_init(jj,k); 
         dggf4_init(jj,k); 
         ggtf1_init(jj,k); 
         dggtf1_init(jj,k); 
         ggtf2_init(jj,k); 
         dggtf2_init(jj,k); 
         ggtf3_init(jj,k);  
         dggtf3_init(jj,k); 
         ggtf4_init(jj,k); 
         dggtf4_init(jj,k)]; 

% =========================================================================



