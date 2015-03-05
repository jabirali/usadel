function [Sxint,Spxint] = deval(sol,xint,idx)
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2004 The MathWorks, Inc.
%   BVP6C Modification
%    Nick Hale  Imperial College London
%    $Date: 12/06/2006 $

if ~isa(sol,'struct')  
  % Try  DEVAL(XINT,SOL)
  temp = sol;
  sol  = xint;
  xint = temp;
end

try
  t = sol.x;    
  y = sol.y;    
catch
  error('MATLAB:deval:SolNotFromDiffEqSolver',...
        '%s must be a structure returned by a differential equation solver.',...
        inputname(1));
end

if nargin < 3
  idx = 1:size(y,1);  % return all solution components
else 
  if any(idx < 0) || any(idx > size(y,1))
    error('MATLAB:deval:IDXInvalidSolComp',...
          'Incorrect solution component requested in %s.',inputname(3));
  end  
end  
idx = idx(:);
  
if isfield(sol,'solver')
  solver = sol.solver;
else
  msg = sprintf('Missing ''solver'' field in the structure %s.',inputname(1));
  if isfield(sol,'yp')
    warning('MATLAB:deval:MissingSolverField',['%s\n         DEVAL will ' ...
            'treat %s as an output of the MATLAB R12 version of BVP4C.'], ...
            msg, inputname(1));  
    solver = 'bvp4c';
  else
    error('MATLAB:deval:NoSolverInStruct', msg);
  end
end

if ~strcmp(solver,'bvp6c')
    if ~ismember(solver,{'ode113','ode15i','ode15s','ode23','ode23s','ode23t',...
                     'ode23tb','ode45','bvp4c','dde23','bvp6c'})
        error('MATLAB:deval:InvalidSolver',...
            'Unrecognized solver name ''%s'' in %s.solver.',solver,inputname(1));
    else 
        % if bvp6c is not the sovler, use existing deval (which ISNT builtin)
        % this is messy, but I cant get it to work otherwise ...
        mydir = cd;
        cd(fullfile(matlabroot, 'toolbox/matlab/funfun/'))
        [Sxint,Spxint] = deval(sol,xint,idx);
        cd(mydir);
        return
    end
end

% Determine the dominant data type
dataType = superiorfloat(sol.x,xint);

Spxint_requested = (nargout > 1);   % Evaluate the first derivative?

n = length(idx);
Nxint = length(xint);
Sxint = zeros(n,Nxint,dataType);
if Spxint_requested
  Spxint = zeros(n,Nxint,dataType);
end

% Make tint a row vector and if necessary, 
% sort it to match the order of t.
tint = xint(:).';  
tdir = sign(t(end) - t(1));
had2sort = any(tdir*diff(tint) < 0);
if had2sort
  [tint,tint_order] = sort(tdir*tint);
  tint = tdir*tint;
end  

% Using the sorted version of tint, test for illegal values.
if (tdir*(tint(1) - t(1)) < 0) || (tdir*(tint(end) - t(end)) > 0)
  error('MATLAB:deval:SolOutsideInterval',...
        ['Attempting to evaluate the solution outside the interval\n'...
         '[%e, %e] where it is defined.\n'],t(1),t(end));
end

% Select appropriate interpolating function.
if ~isfield(sol,'yp') || ~isfield(sol,'ypmid') 
    error('MATLAB:deval:nsufficientData',...
    'bvp6c deval requires the option BVPSET(''SLOPEOUT'',''ON'')');
end
interpfcn = @ntrp6h;

evaluated = 0;
bottom = 1;
while evaluated < Nxint
  
  % Find right-open subinterval [t(bottom), t(bottom+1))
  % containing the next entry of tint. 
  Index = find( tdir*(tint(evaluated+1) - t(bottom:end)) >= 0 );
  bottom = bottom - 1 + Index(end);

  % Is it [t(end), t(end)]?
  at_tend = (t(bottom) == t(end));
      
  if at_tend
    % Use (t(bottom-1) t(bottom)] to interpolate y(t(end)) and yp(t(end)).
    index = find(tint(evaluated+1:end) == t(bottom));
    bottom = bottom - 1;    
  else    
    % Interpolate inside [t(bottom) t(bottom+1)).
    index = find( tdir*(tint(evaluated+1:end) - t(bottom+1)) < 0 );  
  end

  interpdata = { sol.yp(:,bottom) sol.yp(:,bottom+1) ...
        sol.ypmid(:,bottom,:)};               
  
  % Evaluate the interpolant at all points from [t(bottom), t(bottom+1)).
  if Spxint_requested
    [yint,ypint] = feval(interpfcn,tint(evaluated+index),t(bottom),y(:,bottom),...
                         t(bottom+1),y(:,bottom+1),interpdata{:});    
  else  
    yint = feval(interpfcn,tint(evaluated+index),t(bottom),y(:,bottom),...
                 t(bottom+1),y(:,bottom+1),interpdata{:});    
  end

  if at_tend
    bottom = bottom + 1;
  end
  
  % Purify the solution at t(bottom).
  index1 = find(tint(evaluated+index) == t(bottom));
  if ~isempty(index1)    
    yint(idx,index1) = repmat(y(idx,bottom),1,length(index1)); 
  end
  
  % Accumulate the output.
  Sxint(:,evaluated+index) = yint(idx,:);  
  if Spxint_requested
    Spxint(:,evaluated+index) = ypint(idx,:);
  end  
  evaluated = evaluated + length(index);  
end

% For multipoint BVPs, check if solution requested at interface points
multipointBVP = (isequal(solver,'bvp4c') || isequal(solver,'bvp6c')) && any(diff(t)==0);
if multipointBVP
  idxDiscontPoints = find(diff(t)==0);
  for i = idxDiscontPoints
    % Check whether solution requested at the interface
    idxIntPoints = find(tint == t(i));
    n = length(idxIntPoints);
    if n > 0
      % Average the solution across the interface
      nl = '\n         ';
      warning('MATLAB:deval:NonuniqueSolution',...
              ['At the interface XC = %g the solution might not be continuous.'...
               nl 'DEVAL returns the average of the limits from the left and'...
               nl 'from the right of the interface. To approximate the limit'... 
               nl 'values, call DEVAL for XC-EPS(XC) or XC+EPS(XC).\n'],t(i));
      
      Sxint(:,idxIntPoints) = repmat((sol.y(idx,i)+sol.y(idx,i+1))/2,1,n);
      if Spxint_requested
        Spxint(:,idxIntPoints) = repmat((sol.yp(idx,i)+sol.yp(idx,i+1))/2,1,n);
      end  
    end    
  end  
end
  
if had2sort     % Restore the order of tint in the output.
  Sxint(:,tint_order) = Sxint;  
  if Spxint_requested
    Spxint(:,tint_order) = Spxint;  
  end  
end

