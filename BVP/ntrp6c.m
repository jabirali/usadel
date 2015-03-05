function Yint=ntrp6c(f,Xint,x,y,yp,Fmid,varargin)
%NTRP6C  New interpolation helper function for BVP6C.
%   YINT = NTRP6C(F,XINT,SOL) interpolates the bvp6c solution SOL
%   of ode system F to give solution values at new points XINT.
%   Note, it is not necessary to pass F if the slope values SOL.YP
%   and SOL.YPMID are stored in the structure SOL.
%
%   The function may be called without using bvp6c SOL structure using
%   YINT = NTRP6C(F,XINT,X,Y) where X and Y are some obtained solution.
%   Note that if the ODE system contains parameters if F=F(X,Y,P1,P2...)
%   then user should call YINT = NTRP6C(F,XINT,X,Y,[],[],P1,P2,..).
%
%   Note, this function is designed specifically for use via DEVAL.
%
%   See also BVP6C, DEVAL, NTRP6H.

%    Nick Hale  Imperial College London
%    $Date: 12/06/2006 $
expars=[];
if nargin>6     expars=varargin; end
if nargin<6     Fmid=[];         end
if nargin<5     yp=[]; Fmid=[];  end
if nargin<4
    if isstruct(x)
        y=x.y;
        if isfield(x,'yp'),  yp=x.yp;      
        else,                yp=[];           end
        if isfield(x,'ypmid')Fmid=x.ypmid; 
        else,                Fmid=[];         end
        if isfield(x,'parameters') 
            expars=[x.parameters]; 
        end
        if isfield(x,'knownpars')  
            expars=[expars x.knownpars]; 
        end;
        x=x.x;
    else
        error('MATLAB:ntrp6c:NotEnoughInputs',...
            'Insufficient input arguements\n');
    end
end
if (isempty(yp)||isempty(Fmid))&&~isa(f, 'function_handle')
    error('MATLAB:ntrp6c:FNotFunctionHandle',...
        'F must be a function if slope data is not specified.');
end
expars=num2cell(expars);
XL=x(1);
XR=x(end);
N=length(Xint);

if all(XL<=Xint) + all(Xint<=XR)<2    %check that XL<X<XR
    error('MATLAB:ntrp6c:XIntOutsideRange',...
        'X outside interval range [%f,%f]\n',XL,XR);
end

% if all(sort(X)~=X) error('Please order the X[i]'); end

strt=1;
if Xint(1)==XL                   %check if Xint[1]==XL
    I(1)=0; w(1)=0; strt=2;
end
for i=strt:N
    I(i)=0;
    if Xint(i)==XL                  %error if Xint[i]==XL for i!=1
        error('MATLAB:ntrp6c:XLConflict',...
            'Xint[i]=XL for i!=1. Please order the Xint[i]');
    end
    while Xint(i)-x(I(i)+1)>0      %find interval containing each X[i]
        I(i)=I(i)+1;
    end
    w(i)=(Xint(i)-x(I(i)))./(x(I(i)+1)-x(I(i)));
end

J=1;K=1;
if I(1)==0                      %if X[1]=XL
    Yint(:,1)=y(:,1);
    J=2;K=2;
end
while J<=N
    while I(K)-I(J)==0          %find X's in same interval
        K=K+1;
        if K==N+1 
            break;
        end
    end
    if w(K-1)==1                %if X[J]=x[j] then use exact y(:,j)
        Yint(:,(K-1))=y(:,I(J)+1);
        if J<=K-2 
            Yint(:,J:(K-2))=interp(f,x,y,yp,Fmid,expars,w(J:(K-2))',I(J));
        end
    else
        Yint(:,J:(K-1))=interp(f,x,y,yp,Fmid,expars,w(J:(K-1))',I(J));
    end
    J=K;
end

%----------------------------Interpolate----------------------------------------    
function Y=interp(f,x,y,yp,Fmid,expars,w,int)
XL=x(int);  XR=x(int+1);
H=XR-XL;
YL=y(:,int);YR=y(:,int+1);
if isempty(yp)                  %evaluate function f at mesh points
    FL=f(XL,YL,expars{:});      
    FR=f(XR,YR,expars{:});
else                            %mesh slopes are known
    FL=yp(:,int);
    FR=yp(:,int+1);
end
if isempty(Fmid)                %evaluate function f at internal points
    y025 = (54*YL + 10*YR + (9*FL - 3*FR)*H)/64;
    y075 = (10*YL + 54*YR + (3*FL - 9*FR)*H)/64;
    F025 = f(0.25*(3*XL + XR),y025,expars{:}); 
    F075 = f(0.25*(XL + 3*XR),y075,expars{:});
    y05 = 0.5*(YL + YR) - H*(FR-FL+4*(F075-F025))/24;
    F05 = f(0.5*(XL + XR),y05,expars{:});
else                            %interior slopes are known
    F025 = Fmid(:,int,1); 
    F05 =  Fmid(:,int,2); 
    F075 = Fmid(:,int,3); 
end
 
% interpolate at w  
for i=1:length(w)
    wi=w(i);
    Y(:,i)=A66(wi)*YR + A66(1-wi)*YL+     ...
       (XR-XL)*( B66(wi)*FR - B66(1-wi)*FL + ...
      C66(wi)*(F075-F025) + D66(wi)*F05 );
end

function coeff = A66(w)
coeff=w.^2.*polyval([-24 60 -50 15],w);     % w^2*(15-50*w+60*w^2-24*w^3);
function coeff = B66(w)
coeff=w.^2.*polyval([12 -26 19 -5]/3,w);    % w^2/3*(w-1)*(12*w^2-14*w+5);
function coeff = C66(w)
coeff=w.^2.*polyval([-8 16 -8]/3,w);        % -w^2*8/3*(1-w)^2;
function coeff = D66(w)
coeff=w.^2.*polyval([16 -40 32 -8],w);      % w^2*8*(1-w)^2*(2*w-1);