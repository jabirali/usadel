% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Last updated 2015-02-15

classdef SpinVector
    % Define a data structure to describe objects with both a 3x1 vector
    % structure in geometric space and a 2x2 matrix structure in spin
    % space, such as the Pauli vector, and in   general SU(2) vector fields.
    
    % Define internal variables for the data structure
    properties (GetAccess=public, SetAccess=public)
        x = zeros(2);
        y = zeros(2);
        z = zeros(2);
    end
    
    methods
        % Define a constructor for the data structure
        function self = SpinVector(Ax, Ay, Az)
            % Constructs a spin vector from three spin matrices Ax, Ay, Az
            
            % If the arguments are scalars, multiply by identity matrices;
            % If not, just use the input arguments as the components of A.
            if isscalar(Ax)
                self.x = Ax .* eye(2);
            else
                self.x = Ax;
            end
            
            if isscalar(Ay)
                self.y = Ay .* eye(2);
            else
                self.y = Ay;
            end
            
            if isscalar(Az)
                self.z = Az .* eye(2);
            else
                self.z = Az;
            end
        end
        
        
        
        % Overloading of unary operators and functions
        function result = conj(self)
            % This overloads the complex conjugation function
            result = SpinVector(conj(self.x), conj(self.y), conj(self.z));
        end
        
        function result = transpose(self)
            % This overloads the matrix transposition function
            result = SpinVector(self.x.', self.y.', self.z.');
        end
        
        function result = ctranspose(self)
            % This overloads the complex transposition function
            result = SpinVector(self.x', self.y', self.z');
        end
        
        function result = uplus(self)
            % This overloads the unary plus operator
            result = self;
        end
        
        function result = uminus(self)
            % This overloads the unary minus operator
            result = SpinVector(-self.x,-self.y,-self.z);
        end

        
        
        % Overloading of binary operators
        function result = plus(lhs, rhs)
            % This overloads the plus operator for spin vectors
            if isa(lhs, 'SpinVector') && isa(rhs, 'SpinVector')
                result = SpinVector(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
            else
                error('SpinVector:plus', 'You can only add a spin vector to another spin vector.');
            end
        end
        
        function result = minus(lhs, rhs)
            % This overloads the plus operator for spin vectors
            if isa(lhs, 'SpinVector') && isa(rhs, 'SpinVector')
                result = SpinVector(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z);
            else
                error('SpinVector:minus', 'You can only subtract spin vectors from each other.');
            end
        end
        
        function result = times(lhs, rhs)
            % This overloads the arraywise multiplication operator for spin vectors
            if isa(lhs, 'SpinVector') && isa(rhs, 'SpinVector')
                result = SpinVector(lhs.x*rhs.x,lhs.y*rhs.y,lhs.z*rhs.z);
            else
                error('SpinVector:times', 'You can only arraywise multiply two spin vectors.');
            end
        end
        
        function result = mtimes(lhs, rhs)
            % This overloads the matrix multiplication operator for spin vectors
            if isa(lhs, 'SpinVector')
                % Left-hand side is a spin vector
                if isa(rhs, 'SpinVector')
                    % Right-hand side is also a spin vector
                    result = lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
                elseif isscalar(rhs)
                    % Right-hand side is a scalar
                    result = SpinVector(lhs.x*rhs, lhs.y*rhs, lhs.z*rhs);
                elseif isvector(rhs)
                    % Right-hand side is a vector
                    result = lhs.x*rhs(1) + lhs.y*rhs(2) + lhs.z*rhs(3);
                elseif ismatrix(rhs)
                    % Right-hand side is a matrix
                    result = SpinVector(lhs.x*rhs, lhs.y*rhs, lhs.z*rhs);
                else
                    error('SpinVector:mtimes', 'Right-hand argument has to be a spin vector, matrix, vector, or scalar.');
                end
            else
                % Right-hand side is a spin vector
                if isscalar(lhs)
                    % Left-hand side is a scalar
                    result = SpinVector(rhs.x*lhs, rhs.y*lhs, rhs.z*lhs);
                elseif isvector(lhs)
                    % Left-hand side is a vector
                    result = rhs.x*lhs(1) + rhs.y*lhs(2) + rhs.z*lhs(3);
                elseif ismatrix(lhs)
                    % Left-hand side is a matrix
                    result = SpinVector(rhs.x*lhs, rhs.y*lhs, rhs.z*lhs);
                else
                    error('SpinVector:mtimes', 'Left-hand argument has to be a spin vector, matrix, vector, or scalar.');
                end
            end
        end
        
        function result = mpower(lhs,rhs)
            % This overloads the matrix power operator
            if rhs ~= 2
                error('SpinVector:mpower', 'The matrix power of a spin vector has only been defined for an exponent of two.');
            end
            
            result = lhs.x^2 + lhs.y^2 + lhs.z^2;
        end
    end
end
