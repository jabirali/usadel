function sol = bvp6c(ode, bc, solinit, options, varargin) 
%BVP6C  Solve boundary value problems for ODEs by collocation.
%       (6th order extension of BVP4C by Kierzenka and Shampine)
%   SOL = BVP6C(ODEFUN,BCFUN,SOLINIT) integrates a system of ordinary
%   differential equations of the form y' = f(x,y) on the interval [a,b],
%   subject to general two-point boundary conditions of the form
%   bc(y(a),y(b)) = 0. ODEFUN is a function of two arguments: a scalar X
%   and a vector Y. ODEFUN(X,Y) must return a column vector representing
%   f(x,y). BCFUN is a function of two vector arguments. BCFUN(YA,YB) must
%   return a column vector representing bc(y(a),y(b)). SOLINIT is a structure
%   with fields named   
%       x -- ordered nodes of the initial mesh with 
%            SOLINIT.x(1) = a, SOLINIT.x(end) = b
%       y -- initial guess for the solution with SOLINIT.y(:,i)
%            a guess for y(x(i)), the solution at the node SOLINIT.x(i)       
%
%   BVP6C produces a solution that is continuous on [a,b] and has a
%   continuous first derivative there. The solution is evaluated at points
%   XINT using the output SOL of BVP6C and the function DEVAL:
%   YINT = DEVAL(SOL,XINT). The output SOL is a structure with 
%       SOL.x  -- mesh selected by BVP6C
%       SOL.y  -- approximation to y(x) at the mesh points of SOL.x
%       SOL.solver -- 'bvp6c'
%   If specified in BVPSET, SOL may also contain
%       SOL.yp -- approximation to y'(x) at the mesh points of SOL.x
%       sol.ypmid approximations to y'(x) at interior points of SOL.x
%       
%   SOL = BVP6C(ODEFUN,BCFUN,SOLINIT,OPTIONS) solves as above with default
%   parameters replaced by values in OPTIONS, a structure created with the
%   BVPSET function. To reduce the run time greatly, use OPTIONS to supply 
%   a function for evaluating the Jacobian and/or vectorize ODEFUN. 
%   See BVPSET for details and SHOCKBVP for an example that does both.
%
%   SOL = BVP6C(ODEFUN,BCFUN,SOLINIT,OPTIONS,P1,P2...) passes constant, known
%   parameters P1, P2... to the functions ODEFUN and BCFUN, and to all 
%   functions specified in OPTIONS. Use OPTIONS = [] as a place holder if
%   no options are set.   
%   
%   Some boundary value problems involve a vector of unknown parameters p
%   that must be computed along with y(x):
%       y' = f(x,y,p)
%       0  = bc(y(a),y(b),p) 
%   For such problems the field SOLINIT.parameters is used to provide a guess
%   for the unknown parameters. On output the parameters found are returned
%   in the field SOL.parameters. The solution SOL of a problem with one set 
%   of parameter values can be used as SOLINIT for another set. Difficult BVPs 
%   may be solved by continuation: start with parameter values for which you can 
%   get a solution, and use it as a guess for the solution of a problem with 
%   parameters closer to the ones you want. Repeat until you solve the BVP 
%   for the parameters you want.
%
%   The function BVPINIT forms the guess structure in the most common 
%   situations:  SOLINIT = BVPINIT(X,YINIT) forms the guess for an initial mesh X
%   as described for SOLINIT.x and YINIT either a constant vector guess for the
%   solution or a function that evaluates the guess for the solution
%   at any point in [a,b]. If there are unknown parameters, 
%   SOLINIT = BVPINIT(X,YINIT,PARAMS) forms the guess with the vector PARAMS of 
%   guesses for the unknown parameters.  
%
%   BVP6C solves a class of singular BVPs, including problems with 
%   unknown parameters p, of the form
%       y' = S*y/x + f(x,y,p)
%       0  = bc(y(0),y(b),p) 
%   The interval is required to be [0, b] with b > 0. 
%   Often such problems arise when computing a smooth solution of 
%   ODEs that result from PDEs because of cylindrical or spherical 
%   symmetry. For singular problems the (constant) matrix S is
%   specified as the value of the 'SingularTerm' option of BVPSET,
%   and ODEFUN evaluates only f(x,y,p). The boundary conditions
%   must be consistent with the necessary condition S*y(0) = 0 and
%   the initial guess should satisfy this condition.   
%
%   BVP6C can solve multipoint boundary value problems.  For such problems
%   there are boundary conditions at points in [a,b]. Generally these points
%   represent interfaces and provide a natural division of [a,b] into regions.
%   BVP6C enumerates the regions from left to right (from a to b), with indices 
%   starting from 1.  In region k, BVP6C evaluates the derivative as 
%   YP = ODEFUN(X,Y,K).  In the boundary conditions function, 
%   BCFUN(YLEFT,YRIGHT), YLEFT(:,k) is the solution at the 'left' boundary
%   of region k and similarly for YRIGHT(:,k).  When an initial guess is
%   created with BVPINIT(XINIT,YINIT), XINIT must have double entries for 
%   each interface point. If YINIT is a function, BVPINIT calls Y = YINIT(X,K) 
%   to get an initial guess for the solution at X in region k. In the solution
%   structure SOL returned by BVP6C, SOL.x has double entries for each interface 
%   point. The corresponding columns of SOL.y contain the 'left' and 'right' 
%   solution at the interface, respectively. See THREEBVP for an example of
%   solving a three-point BVP.
%
%   Example
%         solinit = bvpinit([0 1 2 3 4],[1 0]);
%         sol = bvp6c(@twoode,@twobc,solinit);
%     solve a BVP on the interval [0,4] with differential equations and 
%     boundary conditions computed by functions twoode and twobc, respectively.
%     This example uses [0 1 2 3 4] as an initial mesh, and [1 0] as an initial 
%     approximation of the solution components at the mesh points.
%         xint = linspace(0,4);
%         yint = deval(sol,xint);
%     evaluate the solution at 100 equally spaced points in [0 4]. The first
%     component of the solution is then plotted with 
%         plot(xint,yint(1,:));
%   For more examples see TWOBVP, FSBVP, SHOCKBVP, MAT4BVP, EMDENBVP, THREEBVP.
%
%   See also BVPSET, BVPGET, BVPINIT, DEVAL, @.
 
%   BVP6C is a finite difference code that implements the Cash-Singhal 6
%   formula. This is a collocation formula and the collocation 
%   polynomial provides a C1-continuous solution that is sixth order
%   accurate uniformly in [a,b]. (For multipoint BVPs, the solution is 
%   C1-continuous within each region, but continuity is not automatically
%   imposed at the interfaces.) Mesh selection and error control are based
%   on the residual of the continuous solution. Analytical condensation is
%   used when the system of algebraic equations is formed.

%   BVP4C
%    Jacek Kierzenka and Lawrence F. Shampine
%    Copyright 1984-2003 The MathWorks, Inc. 
%    $Revision: 1.21.4.8 $  $Date: 2003/12/26 18:08:47 $
%   BVP6C Modification
%    Nick Hale  Imperial College London
%    $Date: 12/06/2006 $


%bvp6c extras
if nargin < 4 
    options= [];
end
MaxNewPts = bvpget(options,'MaxNewPts',2);
solnpts = bvpget(options,'Xint',[]);
if isempty(solnpts) 
    slopeout = strcmp(bvpget(options,'SlopeOut','on'),'on');
else
    slopeout = 0;
    if strcmp(bvpget(options,'SlopeOut'),'on');
          warning('MATLAB:bvp6c:SolnPtsNoSlopeOut', ...
              ['SOL will not contain slope data \n', ...
              '(solution points were specified in BVPSET)']);
    end
end

[m,n]=size(solnpts);
if min(m,n) > 1
    error('MATLAB:bvp6c:SolnPtsNotVector',...
        'The required solution points must be a vector (or []).');
end
if m>n 
    solnpts=solnpts';
end

% check input parameters
if nargin < 3 
  error('MATLAB:bvp6c:NotEnoughInputs', 'Not enough input arguments.')
end

if ~isstruct(solinit)
  error('MATLAB:bvp6c:SolinitNotStruct',...
        'The initial profile must be provided as a structure.')
elseif ~isfield(solinit,'x')
  error('MATLAB:bvp6c:NoXInSolinit',...
        ['The field ''x'' not present in ''' inputname(3) '''.'])
elseif ~isfield(solinit,'y')
  error('MATLAB:bvp6c:NoYInSolinit',...
        ['The field ''y'' not present in ''' inputname(3) '''.'])
end

x = solinit.x(:)';     % row vector
y = solinit.y;   

if isempty(x) || (length(x)<2)
  error('MATLAB:bvp6c:SolinitXNotEnoughPts',...
        ['''' inputname(3) '.x'' must contain at least the two end points.'])
else
  N = length(x);      % number of mesh points
end

xdir = sign(x(end)-x(1));
if any(xdir * diff(x) < 0)
  error ('MATLAB:bvp6c:SolinitXNotMonotonic',...
         ['The entries in ''' inputname(3) ...
          '.x'' must increase or decrease.']);
end

% a multi-point BVP?
mbcidx = find(diff(x) == 0);  % locate internal interfaces
ismbvp = ~isempty(mbcidx);

if isempty(y)
  error('MATLAB:bvp6c:SolinitYEmpty',...
        ['No initial guess provided in ''' inputname(3) '.y''.'])
end
if length(y(1,:)) ~= N
  error('MATLAB:bvp6c:SolXSolYSizeMismatch',...
        ['''' inputname(3) '.y'' not consistent with ''' ...
         inputname(3) '.x''.'])
end

n = size(y,1);	% size of the DE system
nN = n*N;

% stats
nODEeval = 0; 
nBCeval = 0;  

% set the options
if nargin<4
  options = [];
end

% parameters
knownPar = varargin;
unknownPar = isfield(solinit,'parameters');
if unknownPar    
  npar = length(solinit.parameters(:));  
  ExtraArgs = [{solinit.parameters(:)} knownPar];
else
  npar = 0;      
  ExtraArgs = knownPar;    
end               
nExtraArgs = length(ExtraArgs);

% Check the argument functions
if ismbvp
  nregions = length(mbcidx) + 1;
  Lidx = [1, mbcidx+1]; 
  Ridx = [mbcidx, length(x)];
  nBCs = n * nregions + npar;
  testBC = feval(bc,y(:,Lidx),y(:,Ridx),ExtraArgs{:});
else
  nBCs = n + npar;
  testBC = feval(bc,y(:,1),y(:,end),ExtraArgs{:});
end  
nBCeval = nBCeval + 1;
if length(testBC) ~= nBCs
  error('MATLAB:bvp6c:BCfunOuputSize',...
        ['The boundary condition function should return a column vector of' ...
        ' length %i'],nBCs);  
end
if ismbvp 
  region = 1;    % test region 1 only
  testODE = feval(ode,x(1),y(:,1),region,ExtraArgs{:});  
else  
  testODE = feval(ode,x(1),y(:,1),ExtraArgs{:});   
end
nODEeval = nODEeval + 1;
if length(testODE) ~= n
  error('MATLAB:bvp6c:ODEfunOutputSize',...
        ['The derivative function should return a column vector of' ...
        ' length %i'],n);  
end

% Get options and set the defaults
rtol = bvpget(options,'RelTol',1e-3);      
if (length(rtol) ~= 1) || (rtol<=0)
  error('MATLAB:bvp6c:RelTolNotPos', 'RelTol must be a positive scalar.');  
end
if rtol < 100*eps
  rtol = 100*eps;
  warning('MATLAB:bvp6c:RelTolIncrease','RelTol has been increased to %g.',rtol);
end  
atol = bvpget(options,'AbsTol',1e-6);
if length(atol) == 1
  atol = atol(ones(n,1));
else
  if length(atol) ~= n
    error('MATLAB:bvp6c:SizeAbsTol',['Solving %s requires a scalar AbsTol, ' ...
           'or a vector AbsTol of length %d.'],funstring(ode),n);
  end  
  atol = atol(:);
end
if any(atol<=0)
  error('MATLAB:bvp6c:AbsTolNotPos', 'AbsTol must be positive.');
end  

threshold = atol/rtol;

% analytical Jacobians
Fjac = bvpget(options,'FJacobian');
BCjac = bvpget(options,'BCJacobian');  
averageJac = isempty(Fjac); 

Nmax = bvpget(options,'Nmax',floor(2500/n));
printstats = strcmp(bvpget(options,'Stats','off'),'on');

% vectorized (with respect to x and y)
xyVectorized = strcmp(bvpget(options,'Vectorized','off'),'on');
if xyVectorized     
  vectVars = [1,2]; % input to odenumjac
else
  vectVars = [];
end

% Deal with a singular BVP.
singularBVP = false;
S = bvpget(options,'SingularTerm',[]);
if ~isempty(S)
  if (x(1) ~= 0) || (x(end) <= x(1))
    error('MATLAB:bvp6c:SingBVPInvalidInterval',...
          'Singular BVP must be on interval [0, b] with 0 < b.')
  end
  singularBVP = true;
  % Compute matrix for imposing necessary BC, Sy(0) = 0,
  % and impose on guess for solution.
  PBC = eye(size(S)) - pinv(S)*S;
  y(:,1) = PBC*y(:,1);
  % Compute matrix for proper definition of y'(0).
  PImS = pinv(eye(size(S)) - S);   
  % Augment ExtraArgs with information about singular BVP.
  ExtraArgs = [ ExtraArgs {ode Fjac S PImS}];
  ode = @Sode;
  if ~isempty(Fjac)
    if npar > 0
      Fjac = @SFjacpar;
    else
      Fjac = @SFjac;
    end
  end
end

maxNewtIter = 4;  
maxProbes = 4;    % weak line search 
needGlobJac = true;

% Adjust the warnings.
warnstat(1) = warning('query','MATLAB:singularMatrix');
warnstat(2) = warning('query','MATLAB:nearlySingularMatrix');
warnoff = warnstat;
warnoff(1).state = 'off';
warnoff(2).state = 'off';

done = false;  
% THE MAIN LOOP:
while ~done    
  if unknownPar
    Y = [y(:);ExtraArgs{1}];
  else
    Y =  y(:);
  end
    
  [RHS,Xmid,Ymid,yp,Fmid,NF] = colloc_RHS(n,x,Y,ode,bc,npar,xyVectorized, ...
      mbcidx,nExtraArgs,ExtraArgs); 	
  nODEeval = nODEeval + NF;
  nBCeval = nBCeval + 1;

  for iter=1:maxNewtIter
    if needGlobJac
      % setup and factor the global Jacobian            
      [dPHIdy,NF,NBC] = colloc_Jac(n,x,Xmid,Y,Ymid,yp,Fmid,ode,bc,Fjac,BCjac,npar,...
                                   vectVars,averageJac,mbcidx,nExtraArgs,ExtraArgs); 
      needGlobJac = false;                     
      nODEeval = nODEeval + NF;                
      nBCeval = nBCeval + NBC;      
      % explicit row scaling        
      wt = max(abs(dPHIdy),[],2);      
      if any(wt == 0) || ~all(isfinite(nonzeros(dPHIdy)))
        singJac = true;
      else
        scalMatrix = spdiags(1 ./ wt,0,nN+npar,nN+npar);
        dPHIdy = scalMatrix * dPHIdy; 
        [L,U,P] = lu(dPHIdy);        
        singJac = check_singular(dPHIdy,L,U,P,warnstat,warnoff);
      end
      if singJac
        msg = 'Unable to solve the collocation equations -- a singular Jacobian encountered';
        error('MATLAB:bvp6c:SingJac', msg);    
      end 
      scalMatrix = P * scalMatrix;                             
    end
                    
    % find the Newton direction    
    delY = U\(L\( scalMatrix * RHS ));    
    distY = norm(delY);        

    % weak line search with an affine stopping criterion
    lambda = 1;
    for probe = 1:maxProbes     
      Ynew = Y - lambda*delY;                

      if singularBVP     % Impose necessary BC, Sy(0) = 0.
        Ynew(1:n) = PBC*Ynew(1:n);
      end

      if unknownPar  
        ExtraArgs{1} = Ynew(nN+1:end);
      end
      
      [RHS,Xmid,Ymid,yp,Fmid,NF] = colloc_RHS(n,x,Ynew,ode,bc,npar,xyVectorized,mbcidx,nExtraArgs,ExtraArgs);
      
      nODEeval = nODEeval + NF;
      nBCeval = nBCeval + 1;

      distYnew = norm(U \ (L \ (scalMatrix * RHS)));

      if (distYnew < 0.9*distY) 
        break		
      else
        lambda = 0.5*lambda;		
      end
    end  
    
    needGlobJac = (distYnew > 0.1*distY);
          
    if distYnew < 0.1*rtol   
      break
    else
      Y = Ynew;
    end   
  end  
  
 y = reshape(Ynew(1:nN),n,N); % yp, ExtraArgs, and RHS are consistent with y
  
 [res,NF,Yip05,Fmid(:,:,2)] = residual(ode,x,y,yp,Fmid,RHS,threshold,xyVectorized,nBCs, ...
                      mbcidx,ExtraArgs);
  nODEeval = nODEeval + NF; 

  if max(res) < rtol 
    done = true;
  else 	
    % redistribute the mesh
    [N,x,y,mbcidx] = new_profile(n,x,y,Yip05,yp,Fmid,res,mbcidx,rtol,Nmax,MaxNewPts);
    if N > Nmax    
      warning('MATLAB:bvp6c:RelTolNotMet', ...
          [ 'Unable to meet the tolerance without using more than %d '...
            'mesh points. \n The last mesh of %d points and ' ...
            'the solution are available in the output argument. \n ', ...
            'The maximum residual is %g, while requested accuracy is %g.'], ...
            Nmax,length(x),max(res),rtol); 
      sol.solver = 'bvp6c'; 
      sol.maxres = max(res);
      if ~isempty(solnpts)
          sol.note = 'max residual before interpolation';
          sol.x=solnpts;
          sol.y=ntrp6c(ode,sol.x,x,y,yp,Fmid);
      else
          sol.x = x;
          sol.y = y;
          if slopeout
              sol.yp = yp;
              sol.ypmid = Fmid;
          end
      end
      if unknownPar
        sol.parameters = ExtraArgs{1};
      end  
      if ~isempty(knownPar)
        sol.knownpars = knownPar{1};
      end 
      sol.fevals = nODEeval;
      sol.meshexceeded= 1;
      
      if printstats 
          fprintf('The solution was obtained on a mesh of %g points.\n',N);
          fprintf('The maximum residual is %10.3e. \n',sol.maxres); 
          fprintf('There were %g calls to the ODE function. \n',nODEeval); 
          fprintf('There were %g calls to the BC function. \n',nBCeval); 
      end
      
      return  
    end     
    
    nN = n*N;
  
    needGlobJac = true;	
    
  end  

end    % while

% Output
sol.solver = 'bvp6c'; 
if ~isempty(solnpts)
    sol.note = 'max residual before interpolation';
    sol.x=solnpts;
    sol.y=ntrp6c(ode,sol.x,x,y,yp,Fmid);
%     ,ExtraArgs{1},[knownPar{:}]);
else
    sol.x = x;
    sol.y = y;
    if slopeout
      sol.yp = yp;
      sol.ypmid = Fmid;
    end
end
if unknownPar
  sol.parameters = ExtraArgs{1};
end  
if ~isempty(knownPar)
%     sol.knownpars = [knownPar{:}];
    sol.knownpars = [{knownPar(:)}];
end

sol.meshexceeded = 0;
sol.fevals = nODEeval;

% Stats
if printstats 
  fprintf('The solution was obtained on a mesh of %g points.\n',N);
  fprintf('The maximum residual is %10.3e. \n',max(res)); 
  fprintf('There were %g calls to the ODE function. \n',nODEeval); 
  fprintf('There were %g calls to the BC function. \n',nBCeval); 
end


%--------------------------------------------------------------------------

function f = condaux(flag,X,L,U,P)
%CONDAUX  Auxiliary function for estimation of condition.

switch flag
case 'dim'
    f = max(size(L));
case 'real'
    f = 1;
case 'notransp'
    f = U \ (L \ (P * X));
case 'transp'
    f = P * (L' \ (U' \ X));
end


%--------------------------------------------------------------------------

function singular = check_singular(A,L,U,P,warnstat,warnoff)
%CHECK_SINGULAR  Check A (L*U*P) for 'singularity'; mute certain warnings.

[lastmsg,lastid] = lastwarn('');

warning(warnoff);

Ainv_norm = normest1(@condaux,[],[],L,U,P);

warning(warnstat);

[msg,msgid] = lastwarn;
if isempty(msg) 
  lastwarn(lastmsg,lastid);
else
  singular = true;
  for i = 1:length(warnstat)
    if strcmp(msgid,warnstat(i).identifier)
      return;
    end  
  end      
end

singular = (Ainv_norm*norm(A,inf)*eps > 1e-5);

%---------------------------------------------------------------------------

function [Sx, Spx] = interp_Hermite(w,h,y,yp,yp_ip025,yp_ip05,yp_ip075)
%INTERP_HERMITE  use the 6th order Hermite Interpolant presented by Cash
    %and Wright to find y and y' at abscissas for Quadrature.
N=size(y,2);
diagscal = spdiags(h',0,N-1,N-1);

Sx=   A66(w)*y(:,2:N) + A66(1-w)*y(:,1:N-1)   + ...
    ( B66(w)*yp(:,2:N) - B66(1-w)*yp(:,1:N-1) + ...
      C66(w)*(yp_ip075-yp_ip025) + D66(w)*yp_ip05 )*diagscal; 
       
if nargout >1
    diagscal = spdiags(1./h',0,N-1,N-1); 
    Spx=( Ap66(w)*y(:,2:N)  - Ap66(1-w)*y(:,1:N-1) )*diagscal + ...
        ( Bp66(w)*yp(:,2:N) + Bp66(1-w)*yp(:,1:N-1) + ...
          Cp66(w)*(yp_ip075-yp_ip025) + Dp66(w)*yp_ip05 );
end

function coeff = A66(w)
coeff=w.^2.*polyval([-24 60 -50 15],w);     % w^2*(15-50*w+60*w^2-24*w^3);
function coeff = B66(w)
coeff=w.^2.*polyval([12 -26 19 -5]/3,w);    % w^2/3*(w-1)*(12*w^2-14*w+5);
function coeff = C66(w)
coeff=w.^2.*polyval([-8 16 -8]/3,w);        % -w^2*8/3*(1-w)^2;
function coeff = D66(w)
coeff=w.^2.*polyval([16 -40 32 -8],w);      % w^2*8*(1-w)^2*(2*w-1);

function coeff = Ap66(w)
coeff=w.*polyval([-120 240 -150 30],w);     %w*(30-150*w+240*w^2-120*w^3);
function coeff = Bp66(w)
coeff=w.*polyval([20 -104/3 19 -10/3],w);   %w*(w*(20*w^2+19)-(104*w^2+10)/3);
function coeff = Cp66(w)
coeff=-16/3*w.*polyval([2 -3 1],w);         %-16/3*w*(1-3*w+2*w^2);
function coeff = Dp66(w)
coeff=w.*polyval([80 -160 96 -16],w);       %w*(80*w^3-160*w^2+96*w-16);

%---------------------------------------------------------------------------

function [res,nfcn,Y05,Yp05] = residual(Fcn, x, y, yp, Fmid, RHS, threshold, ...
                                xyVectorized, nBCs, mbcidx, ExtraArgs)
%RESIDUAL  Compute L2-norm of the residual using 7-point Lobatto quadrature.

% multi-point BVP support
ismbvp = ~isempty(mbcidx);
nregions = length(mbcidx) + 1;
Lidx = [1, mbcidx+1];
Ridx = [mbcidx, length(x)];
if ismbvp
  FcnArgs = {0,ExtraArgs{:}};  % pass region idx
else  
  FcnArgs = ExtraArgs;
end

lob(2)= 0.0848880518607165;
lobw(2)=0.276826047361566;
lob(3)= 0.265575603264643;
lobw(3)=0.431745381209863;

lob(5)= 0.734424396735357;
lobw(5)=lobw(3);
lob(6)= 0.9151119481392835;
lobw(6)=lobw(2);

[n,N] = size(y);

Yp05=zeros(n,N-1);    %output more accurate yp_ip05
Y05=zeros(n,N-1);    %output more accurate yp_ip05

res = [];
nfcn = 0;

NewtRes = zeros(n,N-1);
% Do not populate the interface intervals for multi-point BVPs.
intidx = setdiff(1:N-1,mbcidx);
NewtRes(:,intidx) = reshape(RHS(nBCs+1:end),n,[]);

for region = 1:nregions
  if ismbvp
    FcnArgs{1} = region;    % Pass the region index to the ODE function.
  end  
  
  xidx = Lidx(region):Ridx(region);    % mesh point index  
  Nreg = length(xidx);
  xreg = x(xidx); 
  yreg = y(:,xidx);
  ypreg = yp(:,xidx);
  hreg = diff(xreg);
  
  iidx = xidx(1:end-1);   % mesh interval index
  Nint = length(iidx);
  thresh = threshold(:,ones(1,Nint));
  
  yp_ip025=Fmid(:,iidx,1);
  yp_ip075=Fmid(:,iidx,3);
  
  diagscal = spdiags(hreg',0,Nint,Nint);
  xip05 = 0.5*(xreg(iidx) + xreg(iidx+1));
  %more accurate estimate of y_ip05 than used in Cash-Singhal 
  Yip05 = 0.5*(yreg(:,iidx+1)+ yreg(:,iidx)) - ...
      (ypreg(:,iidx+1)-yp(:,iidx)+4*(yp_ip075-yp_ip025))*diagscal/24;
  if xyVectorized 
    Yp_ip05 = feval(Fcn,xip05,Yip05,FcnArgs{:}); 
    nfcn = nfcn + 1;    
  else % not vectorized 
    Yp_ip05 = zeros(n,Nint);
    for i = 1:Nint
      Yp_ip05(:,i) = feval(Fcn,xip05(i),Yip05(:,i),FcnArgs{:}); 
    end  
    nfcn = nfcn + Nint;
  end
  Y05(:,iidx)=Yip05;
  Yp05(:,iidx)=Yp_ip05;
 
  res_reg=zeros(1,Nint); 
 
  %Sum contributions from other points
  for j=[2,3,5,6]
      xLob = xreg(1:Nreg-1) + lob(j)*hreg;
      [yLob,ypLob] = interp_Hermite(lob(j),hreg,yreg,ypreg,yp_ip025,Yp_ip05,yp_ip075);      
      if xyVectorized
        fLob = feval(Fcn,xLob,yLob,FcnArgs{:});
        nfcn = nfcn + 1;
      else
        fLob = zeros(n,Nint);
        for i = 1:Nint
          fLob(:,i) = feval(Fcn,xLob(i),yLob(:,i),FcnArgs{:});
        end 
        nfcn = nfcn + Nint;
      end    
      temp = (ypLob - fLob) ./ max(abs(fLob),thresh);
      res_reg = res_reg + lobw(j)*dot(temp,temp,1);
  end 
  
  % scaling
  res_reg = sqrt( abs(hreg/2) .* res_reg);

  res(iidx) = res_reg;  
end        
% %--------------------------------------------------------------------------

function [NN,xx,yy,mbcidxnew] = new_profile(n,x,y,y05,F,Fmid,res,mbcidx,rtol,Nmax,MaxNewPts)
%NEW_PROFILE  Redistribute mesh points and approximate the solution.

% multi-point BVP support
nregions = length(mbcidx) + 1;
Lidx = [1, mbcidx+1];
Ridx = [mbcidx, length(x)];

mbcidxnew = [];

xx = [];
yy = [];
NN = 0;

for region = 1:nregions

  xidx = Lidx(region):Ridx(region);  % mesh point index  
  xreg = x(xidx);
  yreg = y(:,xidx);
  Freg = F(:,xidx);
  hreg = diff(xreg);
  Nreg = length(xidx);

  iidx = xidx(1:end-1);    % mesh interval index
  resreg = res(iidx); 
  
  F025reg = Fmid(:,iidx,1);
  F05reg  = Fmid(:,iidx,2);
  F075reg = Fmid(:,iidx,3);
  
  i1 = find(resreg > rtol);
  i2 = find(resreg(i1) > 100*rtol);
  NNmax = Nreg + length(i1) + length(i2);
  xxreg = zeros(1,NNmax);
  yyreg = zeros(n,NNmax);
  last_int = Nreg - 1;
  
  xxreg(1) = xreg(1);
  yyreg(:,1) = yreg(:,1);
  NNreg = 1;
  i = 1;
  while i <= last_int          
    if resreg(i) > rtol     % introduce points   
      if resreg(i) > 100*rtol   
        Ni = MaxNewPts;
        hi = hreg(i) / (Ni+1);
        j = 1:Ni;
        xxreg(NNreg+j) = xxreg(NNreg) + j*hi;
        for j=1:Ni
            yyreg(:,NNreg+j) = interp_Hermite(j/(Ni+1),hreg(i),yreg(:,i:i+1),...
                     Freg(:,i:i+1),F025reg(:,i),F05reg(:,i),F075reg(:,i));
        end
      else
        Ni = 1;
        xxreg(NNreg+1) = xxreg(NNreg) + hreg(i)/2;
        yyreg(:,NNreg+1) = y05(:,i);          
      end  
      NNreg = NNreg + Ni;
    else                 % try removing points 
      if (i <= last_int-2) && (max(resreg(i+1:i+2)) < rtol)         
        hnew = (hreg(i)+hreg(i+1)+hreg(i+2))/2;
        C1 = resreg(i)/(hreg(i)/hnew)^(11/2);
        C2 = resreg(i+1)/(hreg(i+1)/hnew)^(11/2);
        C3 = resreg(i+2)/(hreg(i+2)/hnew)^(11/2);
        pred_res = max([C1,C2,C3]); 
        
        if pred_res < 0.05 * rtol   % replace 3 intervals with 2 
          xxreg(NNreg+1) = xxreg(NNreg) + hnew;
          yyreg(:,NNreg+1) = interp_Hermite(0.5,hreg(i),yreg(:,(i+1):(i+2)),...
                     Freg(:,(i+1):(i+2)),F025reg(:,i+1),F05reg(:,i+1),F075reg(:,i+1));
          NNreg = NNreg + 1;  
          i = i + 2;    
        end        
      end      
    end
    NNreg = NNreg + 1;
    xxreg(NNreg) = xreg(i+1);   % preserve the next mesh point
    yyreg(:,NNreg) = yreg(:,i+1);
    i = i + 1;
  end  
   
  NN = NN + NNreg;
  if (NN > Nmax)
    % return the previous solution 
    xx = x; 
    yy = y;   
    mbcidxnew = mbcidx;
    break
  else 
    xx = [xx, xxreg(1:NNreg)];
    yy = [yy, yyreg(:,1:NNreg)];
    if region < nregions    % possible only for multipoint BVPs
      mbcidxnew = [mbcidxnew, NN];
    end
  end
end

%--------------------------------------------------------------------------

function [Phi,Xmid,Ymid,F,Fmid,nfcn] = colloc_RHS(n, x, Y, Fcn, Gbc, npar, ...
                                        xyVectorized, mbcidx, nExtraArgs, ExtraArgs)
%COLLOC_RHS  Evaluate the system of collocation equations Phi(Y).  
%   The derivative approximated at the midpoints and returned in Fmid is
%   used to estimate the residual. 

% multi-point BVP support
ismbvp = ~isempty(mbcidx);
nregions = length(mbcidx) + 1;
Lidx = [1, mbcidx+1];
Ridx = [mbcidx, length(x)];
if ismbvp
  FcnArgs = {0,ExtraArgs{:}};    % Pass the region index to the ODE function.
else  
  FcnArgs = ExtraArgs;
end

y = reshape(Y(1:end-npar),n,[]);  

[n,N] = size(y);
nBCs = n*nregions + npar;

F = zeros(n,N);
Fmid = zeros(n,N-1,3);    % include interface intervals
Phi = zeros(nBCs+n*(N-nregions),1);    % exclude interface intervals

% Boundary conditions
% Do not pass info about singular BVPs in ExtraArgs to BC function.
Phi(1:nBCs) = feval(Gbc,y(:,Lidx),y(:,Ridx),ExtraArgs{1:nExtraArgs});    
phiptr = nBCs;    % active region of Phi

for region = 1:nregions
  if ismbvp
    FcnArgs{1} = region;
  end   
  xidx = Lidx(region):Ridx(region);   % mesh point index
  Nreg = length(xidx);
  xreg = x(xidx); 
  yreg = y(:,xidx);
  
  iidx = xidx(1:end-1);   % mesh interval index
  Nint = length(iidx);
  
  % derivative at the mesh points
  if xyVectorized
    Freg = feval(Fcn,xreg,yreg,FcnArgs{:});  
    nfcn = 1;
  else
    Freg = zeros(n,Nreg);
    for i = 1:Nreg
      Freg(:,i) = feval(Fcn,xreg(i),yreg(:,i),FcnArgs{:});  
    end  
    nfcn = Nreg;
  end
  
  %mesh point data
  h = diff(xreg);
  H = spdiags(h(:),0,Nint,Nint);
  xi = xreg(1:end-1);
  yi = yreg(:,1:end-1);
  xip1 = xreg(2:end);
  yip1 = yreg(:,2:end);
  Fi = Freg(:,1:end-1);
  Fip1 = Freg(:,2:end);

  %interior points & derivative
  xip025 = 0.25*(3*xi + xip1);
  xip075 = 0.25*(xi + 3*xip1);
  yip025 = (54*yi + 10*yip1 + (9*Fi - 3*Fip1)*H)/64;
  yip075 = (10*yi + 54*yip1 + (3*Fi - 9*Fip1)*H)/64;
  if xyVectorized 
    Fip025 = feval(Fcn,xip025,yip025,FcnArgs{:}); 
    Fip075 = feval(Fcn,xip075,yip075,FcnArgs{:});
    nfcn = nfcn + 3;    
  else % not vectorized 
    Fip025 = zeros(n,Nint); 
    Fip075 = zeros(n,Nint);
    for i = 1:Nint
      Fip025(:,i) = feval(Fcn,xip025(i),yip025(:,i),FcnArgs{:});
      Fip075(:,i) = feval(Fcn,xip075(i),yip075(:,i),FcnArgs{:});
    end  
    nfcn = nfcn + 2*Nint;
  end
  
  %mid points & derivative  
  xip05 = 0.5*(xi + xip1);
  yip05 = 0.5*(yi + yip1) - (5*Fi - 16*Fip025 + 16*Fip075 - 5*Fip1)*(H/24);
  if xyVectorized 
    Fip05 = feval(Fcn,xip05,yip05,FcnArgs{:}); 
    nfcn = nfcn + 1;    
  else % not vectorized 
    Fip05 = zeros(n,Nint);
    for i = 1:Nint
      Fip05(:,i) = feval(Fcn,xip05(i),yip05(:,i),FcnArgs{:}); 
    end  
    nfcn = nfcn + Nint;
  end
    
  % the Cash-Singhal formula 
  Phireg = yip1 - yi - (7*Fi + 32*Fip025 + 12*Fip05 + 32*Fip075 + 7*Fip1)*(H/90);
 
  % output assembly 
  Phi(phiptr+1:phiptr+n*Nint) = Phireg(:);
  phiptr = phiptr + n*Nint;
   
  Xmid(:,iidx,1) = xip025;
  Xmid(:,iidx,2) = xip05;
  Xmid(:,iidx,3) = xip075;
  
  Ymid(:,iidx,1) = yip025;
  Ymid(:,iidx,2) = yip05;
  Ymid(:,iidx,3) = yip075;
  
  F(:,xidx) = Freg;
  Fmid(:,iidx,1) = Fip025;
  Fmid(:,iidx,2) = Fip05;
  Fmid(:,iidx,3) = Fip075;
  
end

    
%---------------------------------------------------------------------------

function [dPHIdy,nfcn,nbc] = colloc_Jac(n, x, Xmid, Y, Ymid, F, Fmid, ode, bc, ...
                                        Fjac, BCjac, npar, vectVars, averageJac,...
                                        mbcidx, nExtraArgs, ExtraArgs)
%COLLOC_JAC  Form the global Jacobian of collocation eqns.

% multi-point BVP support
ismbvp = ~isempty(mbcidx);
nregions = length(mbcidx) + 1;
Lidx = [1, mbcidx+1];
Ridx = [mbcidx, length(x)];
if ismbvp
  FcnArgs = {0,ExtraArgs{:}};  % pass region idx
else  
  FcnArgs = ExtraArgs;
end

N = length(x);
nN = n*N;               
In = eye(n);

nfcn = 0;   
nbc = 0;
 
y = reshape(Y(1:nN),n,N);

% BC points
ya = y(:,Lidx);
yb = y(:,Ridx);

if isempty(Fjac)  % use numerical approx
  threshval = 1e-6;
  Joptions.diffvar = 2;  % dF(x,y)/dy
  Joptions.vectvars = vectVars;
  Joptions.thresh = threshval(ones(n,1));
  if npar > 0   % unknown parameters
    if ismbvp
      dPoptions.diffvar = 4;  % dF(x,y,region,p)/dp
    else
      dPoptions.diffvar = 3;  % dF(x,y,p)/dp
    end  
    dPoptions.vectvars = vectVars;
    dPoptions.thresh = threshval(ones(npar,1));
    dPoptions.fac = [];        
  end  
end
 
% Collocation equations
nBCs = n*nregions + npar;

rows = nBCs+1:nBCs+n;   % define the action area
cols = 1:n;             % in the global Jacobian

if npar == 0   % no unknown parameters -----------------------------
  
  dPHIdy = spalloc(nN,nN,2*n*nN);  % sparse storage        
  
  if isempty(BCjac)   % use numerical approx
    [dGdya,dGdyb,nbc] = BCnumjac(bc,ya,yb,n,npar,nExtraArgs,ExtraArgs); 
  elseif isa(BCjac,'cell')     % Constant partial derivatives {dGdya,dGdyb}
    dGdya = BCjac{1};
    dGdyb = BCjac{2};        
  else  % use analytical Jacobian 
    [dGdya,dGdyb] = feval(BCjac,ya,yb,ExtraArgs{1:nExtraArgs});
  end   
  
  % Collocation equations 
  for region = 1:nregions
    
    % Left BC
    dPHIdy(1:nBCs,cols) = dGdya(:,(region-1)*n+(1:n));  
        
    if  ismbvp
      FcnArgs{1} = region; % Pass the region index to the ODE function.
    end  
    
    xidx = Lidx(region):Ridx(region);    % mesh point index  
    xreg = x(xidx); 
    yreg = y(:,xidx);
    Freg = F(:,xidx);
    hreg = diff(xreg);
    
    iidx = xidx(1:end-1);    % mesh interval index
    Nint = length(iidx);

    [X1qtrreg, Xmidreg, X3qtrreg] = midptreg(iidx,Xmid);
    [Y1qtrreg, Ymidreg, Y3qtrreg] = midptreg(iidx,Ymid);
    [F1qtrreg, Fmidreg, F3qtrreg] = midptreg(iidx,Fmid);
    
    % Collocation equations
    if isempty(Fjac)  % use numerical approx
            
      Joptions.fac = [];  
      [Ji,Joptions.fac,ignored,nFcalls] = ...
          odenumjac(ode,{xreg(1),yreg(:,1),FcnArgs{:}},Freg(:,1),Joptions);
      nfcn = nfcn+nFcalls;    
      nrmJi = norm(Ji,1);
      
      for i = 1:Nint
        hi = hreg(i);
        % the right mesh point
        xip1 = xreg(i+1);
        yip1 = yreg(:,i+1);
        Fip1 = Freg(:,i+1);        
        [Jip1,Joptions.fac,ignored,nFcalls] = ...
            odenumjac(ode,{xip1,yip1,FcnArgs{:}},Fip1,Joptions); 
        nfcn = nfcn + nFcalls;
        nrmJip1 = norm(Jip1,1);      
        
        %the interior points
        if averageJac && (norm(Jip1 - Ji,1) <= 0.125*(nrmJi + nrmJip1))
          Jip025 = 0.25*(3*Ji + Jip1);
          Jip05 = 0.5*(Ji + Jip1);
          Jip075 = 0.25*(Ji + 3*Jip1);
        else        
          [xip025, xip05, xip075] = midpti(i,X1qtrreg, Xmidreg, X3qtrreg);
          [yip025, yip05, yip075] = midpti(i,Y1qtrreg, Ymidreg, Y3qtrreg);
          [Fip025, Fip05, Fip075] = midpti(i,F1qtrreg, Fmidreg, F3qtrreg);
 
          [Jip025,Joptions.fac,ignored,nFcalls025] = ...
              odenumjac(ode,{xip025,yip025,FcnArgs{:}},Fip025,Joptions);            
          [Jip05,Joptions.fac,ignored,nFcalls05] = ...
              odenumjac(ode,{xip05,yip05,FcnArgs{:}},Fip05,Joptions);            
          [Jip075,Joptions.fac,ignored,nFcalls075] = ...
              odenumjac(ode,{xip075,yip075,FcnArgs{:}},Fip075,Joptions);            
          nfcn = nfcn + nFcalls025 + nFcalls05 + nFcalls075;     
        end 
        
        Jip05Jip025=Jip05*Jip025;
        Jip05Jip075=Jip05*Jip075;
        % assembly
        dPHIdy(rows,cols)=calc_dPHYdy1(hi,In,Ji,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);
        cols = cols + n;
        dPHIdy(rows,cols)=calc_dPHYdy2(hi,In,Jip1,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);
        rows = rows + n;   % next equation     

        Ji = Jip1;
        nrmJi = nrmJip1;
      end        
      
    elseif isa(Fjac,'numeric') % constant Jacobian
      J = Fjac(:,(region-1)*n+(1:n));            
      J2 = J*J;
      J3 = J2*J;
      for i = 1:Nint
        hhreg = hreg(i)*hreg(i);
        hhhreg = hhreg(i)*hreg(i);
        h2J   = hreg(i)/2*J;
        h10J2 = hhhreg/10*J2; 
        h120J3 = hhhreg/120*J3; 
        dPHIdy(rows,cols) = - In - (h2J+h10J2+h120J3);
        cols = cols + n;
        dPHIdy(rows,cols) =   In - (h2J-h10J2+h120J3);
        rows = rows + n;   % next equation
      end  
      
    else % use analytical Jacobian        
      
      Ji = feval(Fjac,xreg(:,1),yreg(:,1),FcnArgs{:});
      
      for i = 1:Nint
        hi = hreg(i);
        % the right mesh point
        xip1 = xreg(i+1);
        yip1 = yreg(:,i+1);
        Jip1 = feval(Fjac,xip1,yip1,FcnArgs{:});    
        
        %the interior points  
        [xip025, xip05, xip075] = midpti(i,X1qtrreg, Xmidreg, X3qtrreg);
        [yip025, yip05, yip075] = midpti(i,Y1qtrreg, Ymidreg, Y3qtrreg);

        Jip025 = feval(Fjac,xip025,yip025,FcnArgs{:});
        Jip05 = feval(Fjac,xip05,yip05,FcnArgs{:});
        Jip075 = feval(Fjac,xip075,yip075,FcnArgs{:});
        
        Jip05Jip025=Jip05*Jip025;
        Jip05Jip075=Jip05*Jip075;
        % assembly
        dPHIdy(rows,cols)=calc_dPHYdy1(hi,In,Ji,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);
        cols = cols + n;
        dPHIdy(rows,cols)=calc_dPHYdy2(hi,In,Jip1,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);
        rows = rows + n;   % next equation  

        Ji = Jip1;
      end
    end
    
    % Right BC
    dPHIdy(1:nBCs,cols) = dGdyb(:,(region-1)*n+(1:n));    
    
    cols = cols + n;
    
  end % regions

else  % there are unknown parameters --------------
  
  dPHIdy = spalloc(nN+npar,nN+npar,(nN+npar)*(2*n+npar));  % sparse storage        
  last_cols = zeros(nN+npar,npar);   % accumulator
  
  if isempty(BCjac)   % use numerical approx
    [dGdya,dGdyb,nbc,dGdpar] = BCnumjac(bc,ya,yb,n,npar,nExtraArgs,ExtraArgs);
  elseif isa(BCjac,'cell')     % Constant partial derivatives {dGdya,dGdyb}
    dGdya  = BCjac{1};
    dGdyb  = BCjac{2};     
    dGdpar = BCjac{3};    
  else  % use analytical Jacobian 
    [dGdya,dGdyb,dGdpar] = feval(BCjac,ya,yb,ExtraArgs{1:nExtraArgs});
  end   
  last_cols(1:nBCs,:) = dGdpar;

  % Collocation equations 
  for region = 1:nregions

    % Left BC
    dPHIdy(1:nBCs,cols) = dGdya(:,(region-1)*n+(1:n));  
        
    if ismbvp
      FcnArgs{1} = region;
    end  
    
    xidx = Lidx(region):Ridx(region);
    xreg = x(xidx); 
    yreg = y(:,xidx);
    Freg = F(:,xidx);
    hreg = diff(xreg);
    
    iidx = xidx(1:end-1);    % mesh interval index
    Nint = length(iidx);
     
    [X1qtrreg, Xmidreg, X3qtrreg] = midptreg(iidx,Xmid);
    [Y1qtrreg, Ymidreg, Y3qtrreg] = midptreg(iidx,Ymid);
    [F1qtrreg, Fmidreg, F3qtrreg] = midptreg(iidx,Fmid);
    
    % Collocation equations
    if isempty(Fjac)  % use numerical approx
            
      Joptions.fac = [];  
      dPoptions.fac = [];
      [Ji,Joptions.fac,dFdpar_i,dPoptions.fac,nFcalls] = ...
          Fnumjac(ode,{xreg(1),yreg(:,1),FcnArgs{:}},Freg(:,1),...
                  Joptions,dPoptions);
      nfcn = nfcn+nFcalls;
      nrmJi = norm(Ji,1);
      nrmdFdpar_i = norm(dFdpar_i,1);
      
      for i = 1:Nint
        hi = hreg(i);
        % the right mesh point
        xip1 = xreg(i+1);
        yip1 = yreg(:,i+1);
        Fip1 = Freg(:,i+1);
        [Jip1,Joptions.fac,dFdpar_ip1,dPoptions.fac,nFcalls] = ...
            Fnumjac(ode,{xip1,yip1,FcnArgs{:}},Fip1,Joptions,dPoptions);    
        nfcn = nfcn + nFcalls; 
        nrmJip1 = norm(Jip1,1);
        nrmdFdpar_ip1 = norm(dFdpar_ip1,1);
        
        %the interior points        
        if averageJac && (norm(Jip1 - Ji,1) <= 0.125*(nrmJi + nrmJip1)) && ...
              (norm(dFdpar_ip1 - dFdpar_i,1) <= 0.125*(nrmdFdpar_i + nrmdFdpar_ip1))
          Jip025 = 0.25*(3*Ji + Jip1);
          Jip05 = 0.5*(Ji + Jip1);
          Jip075 = 0.25*(Ji + 3*Jip1);
          
          dFdpar_ip025 = 0.25*(3*dFdpar_i + dFdpar_ip1);
          dFdpar_ip05 = 0.5*(dFdpar_i + dFdpar_ip1);
          dFdpar_ip075 = 0.25*(dFdpar_i + 3*dFdpar_ip1);
        else         
          [xip025, xip05, xip075] = midpti(i,X1qtrreg, Xmidreg, X3qtrreg);
          [yip025, yip05, yip075] = midpti(i,Y1qtrreg, Ymidreg, Y3qtrreg);
          [Fip025, Fip05, Fip075] = midpti(i,F1qtrreg, Fmidreg, F3qtrreg);

          [Jip025,Joptions.fac,dFdpar_ip025,dPoptions.fac,nFcalls025] = ...
              Fnumjac(ode,{xip025,yip025,FcnArgs{:}},Fip025,Joptions,dPoptions);
          [Jip05,Joptions.fac,dFdpar_ip05,dPoptions.fac,nFcalls05] = ...
              Fnumjac(ode,{xip05,yip05,FcnArgs{:}},Fip05,Joptions,dPoptions);
          [Jip075,Joptions.fac,dFdpar_ip075,dPoptions.fac,nFcalls075] = ...
              Fnumjac(ode,{xip075,yip075,FcnArgs{:}},Fip075,Joptions,dPoptions);

          nfcn = nfcn + nFcalls025 + nFcalls05 + nFcalls075;
        end  

        Jip05Jip025=Jip05*Jip025;
        Jip05Jip075=Jip05*Jip075;
        % assembly
        dPHIdy(rows,cols)=calc_dPHYdy1(hi,In,Ji,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);
        cols = cols + n;
        dPHIdy(rows,cols)=calc_dPHYdy2(hi,In,Jip1,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);
                                 
        last_cols(rows,:) = - ( ...
            hi/90*( 7*(dFdpar_i+dFdpar_ip1)+32*(dFdpar_ip025+dFdpar_ip075)+12*dFdpar_ip05 )     + ...
            hi*hi/180*( Jip025*(9*dFdpar_i-3*dFdpar_ip1) + Jip075*(3*dFdpar_i-9*dFdpar_ip1) + ...
                        Jip05*( 16*(dFdpar_ip025-dFdpar_ip075)+5*(dFdpar_ip1-dFdpar_i) )  )     + ...
            hi*hi*hi/240*( Jip05Jip025*(3*dFdpar_i-dFdpar_ip1)+Jip05Jip075*(3*dFdpar_ip1-dFdpar_i) )  ...
                               );
        rows = rows + n;   % next equation
        
        Ji = Jip1;
        nrmJi = nrmJip1;
        dFdpar_i = dFdpar_ip1;
        nrmdFdpar_i = nrmdFdpar_ip1;
      end
	    
    elseif isa(Fjac,'cell') % Constant partial derivatives {dFdY,dFdp}

      J = Fjac{1}(:,(region-1)*n+(1:n));
      dFdp = Fjac{2}(:,(region-1)*npar+(1:npar));
      J2 = J*J;  
      J3 = J2*J;
      for i = 1:Nint
        hhreg(i) = hreg(i)*hreg(i);
        hhhreg(i) = hhreg(i)*hreg(i);
        h2J   = hreg(i)/2*J;
        h10J2 = hhhreg(i)/10*J2; 
        h120J3 = hhhreg(i)/120*J3; 
        h3_60J2 = hhhreg(i)/60*J2;   
        dPHIdy(rows,cols) = - In - (h2J+h10J2+h120J3);
        cols = cols + n;
        dPHIdy(rows,cols) =   In - (h2J-h10J2+h120J3);       
        
        last_cols(rows,:) = - hreg(i)*dFdp - h3_60J2*DFdp;
        rows = rows+n;   % next equation 
      end

    else % use analytical Jacobian                           
            
      [Ji,dFdpar_i] = feval(Fjac,xreg(1),yreg(:,1),FcnArgs{:});
      
      for i = 1:Nint
        hi = hreg(i);
        % the right mesh point
        xip1 = xreg(i+1);
        yip1 = yreg(:,i+1);
        [Jip1, dFdpar_ip1] = feval(Fjac,xip1,yip1,FcnArgs{:}); 
        
        %the interior points      
        [xip025, xip05, xip075] = midpti(i,X1qtrreg, Xmidreg, X3qtrreg);
        [yip025, yip05, yip075] = midpti(i,Y1qtrreg, Ymidreg, Y3qtrreg);

        [Jip025, dFdpar_ip025] = feval(Fjac,xip025,yip025,FcnArgs{:}); 
        [Jip05, dFdpar_ip05] = feval(Fjac,xip05,yip05,FcnArgs{:}); 
        [Jip075, dFdpar_ip075] = feval(Fjac,xip075,yip075,FcnArgs{:});         
 
 	Jip05Jip025 = Jip05*Jip025;  
        Jip05Jip075 = Jip05*Jip075;  
        
        % assembly
        dPHIdy(rows,cols)=calc_dPHYdy1(hi,In,Ji,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);
        cols = cols + n;
        dPHIdy(rows,cols)=calc_dPHYdy2(hi,In,Jip1,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075);
                                 
        last_cols(rows,:) = - ( ...
            hi/90*( 7*(dFdpar_i+dFdpar_ip1)+32*(dFdpar_ip025+dFdpar_ip075)+12*dFdpar_ip05 )     + ...
            hi*hi/180*( Jip025*(9*dFdpar_i-3*dFdpar_ip1) + Jip075*(3*dFdpar_i-9*dFdpar_ip1) + ...
                        Jip05*( 16*(dFdpar_ip025-dFdpar_ip075)+5*(dFdpar_ip1-dFdpar_i) )  )     + ...
            hi*hi*hi/240*( Jip05Jip025*(3*dFdpar_i-dFdpar_ip1)+Jip05Jip075*(3*dFdpar_ip1-dFdpar_i) )  ...
                               );
        rows = rows + n;   % next equation
        
        Ji = Jip1;
        dFdpar_i = dFdpar_ip1;
      end
    end
          
    % Right BC
    dPHIdy(1:nBCs,cols) = dGdyb(:,(region-1)*n+(1:n));    

    cols = cols + n;	  
  end
  
  dPHIdy(:,end-npar+1:end) = last_cols;  % accumulated
end
    
%--------------------------------------------------------------------------

function dPHIdy=calc_dPHYdy1(hi,In,Ji,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075)

     dPHIdy   = - In -             ( ...
            hi/90*       ( 7*Ji+27*Jip025+6*Jip05+5*Jip075   ) + ...
            hi*hi/360*   ( 27*Jip05Jip025-5*Jip05Jip075)       + ...
            ( hi*hi/360*(18*Jip025-10*Jip05+6*Jip075) + ...
              hi*hi*hi/240*(3*Jip05Jip025-Jip05Jip075) )*Ji      ...
                                   );
%--------------------------------------------------------------------------

function dPHIdy=calc_dPHYdy2(hi,In,Jip1,Jip025,Jip05,Jip075,Jip05Jip025,Jip05Jip075)

      dPHIdy   = In -               (...
            hi/90*       ( 5*Jip025+6*Jip05+27*Jip075+7*Jip1  ) + ...
            hi*hi/360*   ( 5*Jip05Jip025-27*Jip05Jip075 )       - ...
            ( hi*hi/360*(6*Jip025-10*Jip05+18*Jip075) + ...             
              hi*hi*hi/240*(Jip05Jip025-3*Jip05Jip075)  )*Jip1   ...
                                    ); 
%--------------------------------------------------------------------------

function res = bcaux(Ya,Yb,n,bcfun,varargin)
  ya = reshape(Ya,n,[]);
  yb = reshape(Yb,n,[]);
  res = feval(bcfun,ya,yb,varargin{:});

%---------------------------------------------------------------------------

function [dBCdya,dBCdyb,nCalls,dBCdpar] = BCnumjac(bc,ya,yb,n,npar, ...
                                                  nExtraArgs,ExtraArgs)
%BCNUMJAC  Numerically compute dBC/dya, dBC/dyb, and dBC/dpar, if needed.

% Do not pass info about singular BVPs in ExtraArgs to BC function.
bcArgs = {ya(:),yb(:),n,bc,ExtraArgs{1:nExtraArgs}};  
dBCoptions.thresh = repmat(1e-6,length(ya(:)),1);
dBCoptions.fac = [];
dBCoptions.vectvars = []; % BC functions not vectorized

bcVal  = bcaux(bcArgs{:});
nCalls = 1;

dBCoptions.diffvar = 1;
[dBCdya,ignored,ignored1,nbc] = odenumjac(@bcaux,bcArgs,bcVal,dBCoptions);
nCalls = nCalls + nbc;
dBCoptions.diffvar = 2;
[dBCdyb,ignored,ignored1,nbc] = odenumjac(@bcaux,bcArgs,bcVal,dBCoptions);
nCalls = nCalls + nbc;
if npar > 0
  bcArgs = {ya,yb,ExtraArgs{1:nExtraArgs}};  
  dBCoptions.thresh = repmat(1e-6,npar,1);
  dBCoptions.diffvar = 3;  
  [dBCdpar,ignored,ignored1,nbc] = odenumjac(bc,bcArgs,bcVal,dBCoptions);
  nCalls = nCalls + nbc;
end
  
%---------------------------------------------------------------------------

function [dFdy,dFdy_fac,dFdp,dFdp_fac,nFcalls] = Fnumjac(ode,odeArgs,odeVal,...
                                                    Joptions,dPoptions)
%FNUMJAC  Numerically compute dF/dy and dF/dpar. 

[dFdy,dFdy_fac,ignored,dFdy_nfcn] = odenumjac(ode,odeArgs,odeVal,Joptions);  

[dFdp,dFdp_fac,ignored,dFdp_nfcn] = odenumjac(ode,odeArgs,odeVal,dPoptions);

nFcalls = dFdy_nfcn + dFdp_nfcn;
  
%---------------------------------------------------------------------------

function Fval = Sode(x,y,varargin)
%SODE    Evaluate the derivative function for singular BVPs.
%   Information for singular BVP starts at varargin{end-3}.
odefun = varargin{end-3};
Fval = feval(odefun,x,y,varargin{1:end-4});

% singular point, PImS = varargin{end};
idx = find(x == 0);
if ~isempty(idx)
  Fval(:,idx) = varargin{end}*Fval(:,idx);
end

% regular points, S = varargin{end-1}
idx = find(x ~= 0);
if ~isempty(idx)
  Fval(:,idx) = Fval(:,idx) + (varargin{end-1}*y(:,idx)) * ...
                              spdiags(1./x(idx)',0,length(idx),length(idx));
end

%---------------------------------------------------------------------------

function dFdy = SFjac(x,y,varargin)
%SFJAC   Numerically compute dF/dy for singular BVPs.
%   Information for singular BVP starts at varargin{end-3}.
Fjac = varargin{end-2};
if isa(Fjac,'numeric')
  n = size(Fjac,1);
  cols = 1:n;
  % Is this a multipoint BVP?
  nreg = size(Fjac,2)/n;   % rectangular Jacobian?
  if nreg > 1
    region = varargin{1};   
    cols = n*(region-1)+1 : n*region;
  end
  dFdy = Fjac(:,cols);
else  
  dFdy = feval(Fjac,x,y,varargin{1:end-4});
end
if x > 0
  % S = varargin{end-1}
  dFdy = dFdy + varargin{end-1}/x; 
else
  % PImS = varargin{end};
  dFdy = varargin{end}*dFdy;
end

%---------------------------------------------------------------------------

function [dFdy, dFdpar] = SFjacpar(x,y,varargin)
%SFJACPAR  Numerically compute dF/dy and dF/dpar for singular BVPs.
%   Information for singular BVP starts at varargin{end-3}.
Fjac = varargin{end-2};
if isa(Fjac,'cell')  
  n_y = size(Fjac{1},1);
  % Is this a multipoint BVP?
  nreg = size(Fjac{1},2)/n;   % rectangular Jacobian?
  n_p = size(Fjac{2},2)/nreg; % number of parameters 
  cols_y = 1:n_y;
  cols_p = 1:n_p;
  % Is this a multipoint BVP?
  if nreg > 1
    region = varargin{1};   
    cols_y = n_y*(region-1)+1 : n_y*region;
    cols_p = n_p*(region-1)+1 : n_p*region;
  end
  dFdy   = Fjac{1}(:,cols_y);   
  dFdpar = Fjac{2}(:,cols_p);
else  
  [dFdy, dFdpar] = feval(Fjac,x,y,varargin{1:end-4});
end
if x > 0
  % S = varargin{end-1}
  dFdy = dFdy + varargin{end-1}/x;
else
  % PImS = varargin{end};
  dFdy = varargin{end}*dFdy;
  dFdpar = varargin{end}*dFdpar;
end


%--------------------------------------------------------------------------

function [a025, a05, a075] = midptreg(iidx,mid)
a025 = mid(:,iidx,1);
a05 = mid(:,iidx,2);
a075 = mid(:,iidx,3);

function [a025, a05, a075] = midpti(i, b025, b05, b075)
a025 = b025(:,i);
a05 = b05(:,i);
a075 = b075(:,i);
