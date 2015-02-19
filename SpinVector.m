% Define a data structure to describe objects with both a 3x1 vector
% structure in geometric space and a 2x2 matrix structure in spin
% space, such as the Pauli vector, and in general SU(2) vector fields.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-14
% Updated 2015-02-16
    
classdef SpinVector
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the internal variables for the data structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (GetAccess=public, SetAccess=public)
        x = zeros(2);
        y = zeros(2);
        z = zeros(2);
    end

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the internal methods and operator overloading
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function self = SpinVector(Ax, Ay, Az)
            % Constructs a spin vector from three spin matrices Ax, Ay, Az
            self.x = Ax;
            self.y = Ay;
            self.z = Az;
            
            % If the arguments are scalars, multiply by identity matrices;
            % If not, just use the input arguments as the components of A.
            % if isscalar(Ax)
            %     self.x = Ax .* eye(2);
            % else
            %     self.x = Ax;
            % end
            % 
            % if isscalar(Ay)
            %     self.y = Ay .* eye(2);
            % else
            %     self.y = Ay;
            % end
            % 
            % if isscalar(Az)
            %     self.z = Az .* eye(2);
            % else
            %     self.z = Az;
            % end
        end
        
        
        
        % Overloading of display functions
        function display(self)
            name = inputname(1);
            disp(sprintf(':: %s.x:', name));
            disp([self.x]);
            disp(sprintf('\n:: %s.y:', inputname(1)));
            disp([self.y]);
            disp(sprintf('\n:: %s.z:', inputname(1)));
            disp([self.z]);
        end        
        
        
        
        % Overloading of unary operators and functions
        function self = conj(self)
            % This overloads the complex conjugation function
            self.x = conj(self.x);
            self.y = conj(self.y);
            self.z = conj(self.z);
        end
        
        function self = transpose(self)
            % This overloads the matrix transposition function
            self.x = self.x.';
            self.y = self.y.';
            self.z = self.z.';
        end
        
        function self = ctranspose(self)
            % This overloads the complex transposition function
            self.x = self.x';
            self.y = self.y';
            self.z = self.z';
        end
        
        function self = uplus(self)
            % This overloads the unary plus operator
        end
        
        function self = uminus(self)
            % This overloads the unary minus operator
            self.x = -self.x;
            self.y = -self.y;
            self.z = -self.z;
        end

        
        
        % Overloading of binary operators
        function lhs = plus(lhs, rhs)
            % This overloads the plus operator for spin vectors
            lhs.x = lhs.x+rhs.x;
            lhs.y = lhs.y+rhs.y;
            lhs.z = lhs.z+rhs.z;
            
            % if isa(lhs, 'SpinVector') && isa(rhs, 'SpinVector')
            %     result = SpinVector(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
            % else
            %     error('SpinVector.plus: You can only add a spin vector to another spin vector.');
            % end
        end
        
        function lhs = minus(lhs, rhs)
            % This overloads the plus operator for spin vectors
            lhs.x = lhs.x-rhs.x;
            lhs.y = lhs.y-rhs.y;
            lhs.z = lhs.z-rhs.z;
            
            % if isa(lhs, 'SpinVector') && isa(rhs, 'SpinVector')
            %     result = SpinVector(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z);
            % else
            %     error('SpinVector.minus: You can only subtract spin vectors from each other.');
            % end
        end
        
        function lhs = times(lhs, rhs)
            % This overloads the arraywise multiplication operator for spin vectors
            lhs.x = lhs.x*rhs.x;
            lhs.y = lhs.y*rhs.y;
            lhs.z = lhs.z*rhs.z;

            % if isa(lhs, 'SpinVector') && isa(rhs, 'SpinVector')
            %     result = SpinVector(lhs.x*rhs.x,lhs.y*rhs.y,lhs.z*rhs.z);
            % else
            %     error('SpinVector.times: You can only arraywise multiply two spin vectors.');
            % end
        end
        
        function result = mtimes(lhs, rhs)
            % This overloads the matrix multiplication operator for spin vectors

            if isobject(lhs)
                % Left-hand side is a spin vector
                if isobject(rhs)
                    % Right-hand side is a spin vector
                    result = lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
                elseif length(rhs) == 3
                    % Right-hand side is a vector
                    result = lhs.x*rhs(1) + lhs.y*rhs(2) + lhs.z*rhs(3);
                else
                    % Right-hand side is a matrix or scalar
                    lhs.x = lhs.x*rhs;
                    lhs.y = lhs.y*rhs;
                    lhs.z = lhs.z*rhs;
                    result = lhs;
                end
            else
                % Right-hand side is a spin vector, left-hand side is not
                if length(lhs) == 3
                    % Left-hand side is a vector
                    result = rhs.x*lhs(1) + rhs.y*lhs(2) + rhs.z*lhs(3);
                else
                    % Left-hand side is a matrix or scalar
                    rhs.x = rhs.x*lhs;
                    rhs.y = rhs.y*lhs;
                    rhs.z = rhs.z*lhs;
                    result = rhs;
                end
            end
        end
        
        function result = mpower(lhs,rhs)
            % This overloads the matrix power operator. Note that the
            % current definition is only valid for even exponents!
            
            % if rhs ~= 2
            %     error('SpinVector.mpower: The matrix power of a spin vector has only been defined for an exponent of two.');
            % end
            
            result = (lhs.x^2 + lhs.y^2 + lhs.z^2)^(rhs/2);
        end
        
        function lhs = mrdivide(lhs,rhs)
            % This overloads the matrix division operator            
            lhs.x = lhs.x/rhs;
            lhs.y = lhs.y/rhs;
            lhs.z = lhs.z/rhs;
        end
    end
    
   
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define constant properties (available without object instantiation)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    properties (Constant)
          Pauli = SpinVector([0,1;1,0], [0,-i;i,0], [1,0;0,-1]);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define static methods (available without object instantiation)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static)        
        function result = RashbaDresselhaus(strength, angle)
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

        % Define the Rashba--Dresselhaus SU(2) field
        result = SpinVector( strength*(cos(angle)*SpinVector.Pauli.x   ...
                                     + sin(angle)*SpinVector.Pauli.y), ...
                            -strength*(cos(angle)*SpinVector.Pauli.y   ...
                                     + sin(angle)*SpinVector.Pauli.x), ...
                             0 );
        end
    end
end