% Wr    itten by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Based on a similar program by Dr. Sol Jacobsen
% Created 2015-02-15
% Updated 2015-02-15


function MainSF(BCS)
%
% Input:    
%           BCS     'true':  use bulk solution as the initial guess for g;
%                   'false': use zero matrices as the initial guess for g.
%                   In general, you want to set this to 'true' for strong
%                   proximity effect, and 'false' for weak proximity effect.

% Add subdirectories to the Matlab search path
addpath('BVP/')     % Boundary value problem solver

% Define a grid with positions and energy values
[X,E] = meshgrid(linspace(Emin,Emax), linspace(Xmin,Xmax));


% TODO: PCHIP INTERPOLATION


% PARAMETER CHOICES
% All energy units measured in Delta. The x-coordinate is normalized
% against the length of the F region and thus ranges from 0 to 1. There is
% a superconductor at x=0 and at x=1. 

zeta = 3; % The ratio interface/bulk resistance. Set equal for the moment but will need two separate zetas

xi = 30; % SC coherence length in bulk in nm
dF = 15; % F length in nm
%dF = 30;

% The ratio between SC coherence length and dF squared is equal to the
% Thouless energy divided on the SC gap, as seen from the definition of
% Ethf = D/dF^2 and xi = sqrt(D/Delta). 
Ethf = (xi/dF)^2; %Ethf = ((xi/dF)^2)*Delta;...
delta = 1e-03; % add small imaginary value to energies for numerical 
%stability and to model inelastic scattering

% Define the energies you want to solve the Usadel equation for. In terms
% of making predictions for spectroscopic features such as the density of
% states, the relevant energy range is typically from E=0 to E=2*Delta
% since for E>>Delta, we just see a flat DOS (the normal-state DOS). The
% signatures of superconducting correlations occur for energies of order
% Delta. 
stepE = 0.01;
Evec = [-1.3:stepE:1.3];

% Predefine matrices of zero (makes the code a little quicker). The ggf
% matrices will contain the solution whereas the ggf_init matrices contain
% the initial guesses used to solve the Usadel equation. xvec is a vector
% spanning the range 0 to 1 with 100 steps. xx is the length of this
% vector, thus equal to 100.
xvec = linspace(0,1);
xx = length(linspace(0,1));

% gs (as in plural of g) 
gs = 2

ggf1 = zeros(length(Evec),xx);
ggf2 = zeros(length(Evec),xx);
ggf3 = zeros(length(Evec),xx);
ggf4 = zeros(length(Evec),xx);

ggtf1 = zeros(length(Evec),xx);
ggtf2 = zeros(length(Evec),xx);
ggtf3 = zeros(length(Evec),xx);
ggtf4 = zeros(length(Evec),xx);

dggf1 = zeros(length(Evec),xx);
dggf2 = zeros(length(Evec),xx);
dggf3 = zeros(length(Evec),xx);
dggf4 = zeros(length(Evec),xx);

dggtf1 = zeros(length(Evec),xx);
dggtf2 = zeros(length(Evec),xx);
dggtf3 = zeros(length(Evec),xx);
dggtf4 = zeros(length(Evec),xx);

% Predefine matrices of zero for the initial guesses
ggf1_init = zeros(length(Evec),xx);
ggf2_init = zeros(length(Evec),xx);
ggf3_init = zeros(length(Evec),xx);
ggf4_init = zeros(length(Evec),xx);

ggtf1_init = zeros(length(Evec),xx);
ggtf2_init = zeros(length(Evec),xx);
ggtf3_init = zeros(length(Evec),xx);
ggtf4_init = zeros(length(Evec),xx);

dggf1_init = zeros(length(Evec),xx);
dggf2_init = zeros(length(Evec),xx);
dggf3_init = zeros(length(Evec),xx);
dggf4_init = zeros(length(Evec),xx);

dggtf1_init = zeros(length(Evec),xx);
dggtf2_init = zeros(length(Evec),xx);
dggtf3_init = zeros(length(Evec),xx);
dggtf4_init = zeros(length(Evec),xx);


% Define vectors with the values we want to loop over - note no phasevec
% for SF only. 
hxvec = [6,0,6,0,5,0];
hyvec = [4,4,5,5,6,6];
hzvec = [0,6,0,6,0,5];

Rvec = [0,1,4];
Dvec = [0,1,4];

% Define cell arrays of zeros to store the DoS and triplet matrices
NFarray = cell(length(hzvec),length(Rvec));
triparray1 = cell(length(hzvec),length(Rvec));
triparray2 = cell(length(hzvec),length(Rvec));
tripnoprojarray = cell(length(hzvec),length(Rvec));
singarray = cell(length(hzvec),length(Rvec));
reperpcomparray = cell(length(hzvec),length(Rvec));

for f=1:length(hzvec)
    for h=1:length(Rvec)
        display(f);
        display(h);
        for k=1:1:length(Evec)
                
               hx= hxvec(f);
               hy= hyvec(f);
               hz= hzvec(f);

               A = {Axhat,Ayhat,Azhat}; 
                    
               E=Evec(k) + i*delta;
               display(E);
               jj = k;
    
               % Solve the Usadel equation in the F region
               [ggf1(k,:),ggf2(k,:),ggf3(k,:),ggf4(k,:),ggtf1(k,:),ggtf2(k,:),ggtf3(k,:),ggtf4(k,:),dggf1(k,:),dggf2(k,:),dggf3(k,:),dggf4(k,:),dggtf1(k,:),dggtf2(k,:),dggtf3(k,:),dggtf4(k,:)] = Sol_SolveUsadel_SO_SF_ForAli();

                    if k ~= length(Evec)
                    ggf1_init(k+1,:) = ggf1(k,:);
                    ggf2_init(k+1,:) = ggf2(k,:);
                    ggf3_init(k+1,:) = ggf3(k,:);
                    ggf4_init(k+1,:) = ggf4(k,:);
        
                    ggtf1_init(k+1,:) = ggtf1(k,:);
                    ggtf2_init(k+1,:) = ggtf2(k,:);
                    ggtf3_init(k+1,:) = ggtf3(k,:);
                    ggtf4_init(k+1,:) = ggtf4(k,:);
        
                    dggf1_init(k+1,:) = dggf1(k,:);
                    dggf2_init(k+1,:) = dggf2(k,:);
                    dggf3_init(k+1,:) = dggf3(k,:);
                    dggf4_init(k+1,:) = dggf4(k,:);
        
                    dggtf1_init(k+1,:) = dggtf1(k,:);
                    dggtf2_init(k+1,:) = dggtf2(k,:);
                    dggtf3_init(k+1,:) = dggtf3(k,:);
                    dggtf4_init(k+1,:) = dggtf4(k,:);
                    end  
    
                    for j=1:1:length(xvec)
                        g1 = ggf1(k,j);
                        g2 = ggf2(k,j);
                        g3 = ggf3(k,j);
                        g4 = ggf4(k,j);
        
                        gt1 = ggtf1(k,j);
                        gt2 = ggtf2(k,j);
                        gt3 = ggtf3(k,j);
                        gt4 = ggtf4(k,j);
        
                        gammaf = [g1,g2;g3,g4];
                        gammatf = [gt1,gt2;gt3,gt4];
        
                        N = inv(id-gammaf*gammatf);
                        Nt = inv(id-gammatf*gammaf);
        
                        NF(k,j) = 0.5*trace(N*(id+gammaf*gammatf));
        
                        % Also get triplet components to identify LRT
        
                        tripmat=2*N*gammaf;
                        tripnoproj(k,j)=0.5*(tripmat(1,2)+tripmat(2,1));
                        trip1(k,j)=0.5*(-1i*(tripmat(1,1)+tripmat(2,2)));
                        trip2(k,j)=0.5*(tripmat(1,1)-tripmat(2,2)); %Wrong order...
                        singfn(k,j)=0.5*(tripmat(1,2)-tripmat(2,1));
                        
                        NFarray{f,h}=NF;
                        tripnoprojarray{f,h}=tripnoproj;
                        triparray1{f,h}=trip1;
                        triparray2{f,h}=trip2;
                        singarray{f,h}=singfn;
                        
                        tripvec=[trip2(k,j),trip1(k,j),tripnoproj(k,j)];
                        fieldvec=[hx,hy,hz];
                        reperpcomp(k,j)=norm(real(cross(tripvec,fieldvec)))/norm(fieldvec);
                        
                        reperpcomparray{f,h}=reperpcomp;
                       
                        
                    end %end j loop
        end % end k
    end %end h
end %end f

% For plotting individual lines of DoS at a specific point in the N,
% use syntax like plot(Evec,real(NFarray{f,h}(:,2))), where f sets the h,
% g sets the phase and h picks specific R and D values. The (:,2) plots the
% DOS for all energies at x=2 units away from the edge. Middle plotted by
% (:,50) and so on. 
