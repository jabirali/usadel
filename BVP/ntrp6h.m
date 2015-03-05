function [yint,ypint] = ntrp6h(tint,t,y,tnew,ynew,yp,ypnew,Fmid)
%NTRP6H  Interpolation helper function for BVP6C.
%   YINT = NTRP6H(TINT,T,Y,TNEW,YNEW,YP,YPNEW,FMID) evaluates the
%   Cash-Moore, Cash-Singhal6 based linear interpolant at time TINT.
%   TINT may be a scalar or a row vector.   
%   [YINT,YPINT] = NTRP6H(TINT,T,Y,TNEW,YNEW,YP,YPNEW,FMID) returns
%   also the derivative of the interpolating polynomial. 
%   
%   See also BVP6C, DEVAL, NTRP6C

%    Nick Hale  Imperial College London
%    $Date: 12/06/2006 $

h = tnew - t;
w = (tint - t)/h;
for i=1:length(tint)
    yint(:,i) = A66(w(i))*ynew + A66(1-w(i))*y   + ...
            ( B66(w(i))*ypnew - B66(1-w(i))*yp + ...
              C66(w(i))*(Fmid(:,:,3)-Fmid(:,:,1)) + D66(w(i))*Fmid(:,:,2) )*h;         
    if nargout > 1
        ypout(:,i) =( Ap66(w(i))*ynew  - Ap66(1-w(i))*y )/h + ...
               ( Bp66(w(i))*ypnew + Bp66(1-w(i))*yp + ...
                 Cp66(w(i))*(Fmid(:,:,3)-Fmid(:,:,1)) + Dp66(w(i))*Fmid(:,:,2) );  
    end
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

