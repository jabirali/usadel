% Define a data structure to describe the state of the physical system
% for a given position and energy. This is done by describing the
% Riccati parametrized Green's function 'g', it's tilde conjugate 'gt', 
% and their first derivatives 'dg' and 'dgt' for that configuration.
%
% This class is mainly intended for use with differential equation
% solvers, and therefore provides the method 'vectorize' to pack the
% internal variables in a vector format, and constructor State(...)
% that is able to unpack this vector format.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-14
% Updated 2015-02-16


classdef State
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the internal variables of the data structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    properties (GetAccess=public, SetAccess=public)
        g   = zeros(2);
        dg  = zeros(2);
        gt  = zeros(2);
        dgt = zeros(2);
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the internal methods and overloaded operators
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    methods
        function self = State(varargin)
            % This is the default constructor, which takes as its input either:
            %  (i)   four 2x2 matrices, which should correspond to g, dg, gt, dgt;
            %  (ii)  one 16-element complex vector, which should be produced by 'vectorize' method;
            %  (iii) no arguments, i.e. the empty constructor.
            switch nargin
                case 1
                    % If we get one input, then assume that we got a vector
                    % created by 'vectorize', and reverse the procedure
                    
                    args     = reshape(varargin{1}, 2, 8);
                    self.g   = args(:,1:2);
                    self.dg  = args(:,3:4);
                    self.gt  = args(:,5:6);
                    self.dgt = args(:,7:8);
                case 4
                    % If we get four input arguments, then assume that
                    % these correspond to g, dg, gt, and dgt, respectively.
                    % Multiply with a 2x2 identity matrix in case the
                    % input arguments were scalars.
                    
                    self.g   = varargin{1} * eye(2);
                    self.dg  = varargin{2} * eye(2);
                    self.gt  = varargin{3} * eye(2);
                    self.dgt = varargin{4} * eye(2);
                otherwise
                    % In any other case, we assume that we were called as
                    % an empty constructor, so do nothing to the members.
            end
        end
        
        % Overloading of the display function
        function display(self)
            name = inputname(1);
            fprintf(':: %s.g:            [  Riccati parameter g  ]\n', name);
            disp([self.g]);
            fprintf('\n:: %s.dg:           [   Derivative  dg/dz   ]\n', inputname(1));
            disp([self.dg]);
            fprintf('\n:: %s.gt:           [  Riccati parameter g~ ]\n', inputname(1));
            disp([self.gt]);
            fprintf('\n:: %s.dgt:          [   Derivative  dg~/dz  ]\n', inputname(1));
            disp([self.dgt]);
        end

        % Definition of other useful methods
        function result = vectorize(self)
            % Convert the internal data structure to a vector shape
            result = reshape([self.g self.dg self.gt self.dgt], 1, 16);
                     %[self.vectorize_g  self.vectorize_dg ...
                      %self.vectorize_gt self.vectorize_dgt]';
        end
        
        function result = vectorize_g(self)
            % Convert part of the internal data structure to a vector shape
            result = reshape(self.g,  1, 4);
        end
        
        function result = vectorize_dg(self)
            % Convert part of the internal data structure to a vector shape
            result = reshape(self.dg,  1, 4);
        end
        
        function result = vectorize_gt(self)
            % Convert part of the internal data structure to a vector shape
            result = reshape(self.gt,  1, 4);
        end
        
        function result = vectorize_dgt(self)
            % Convert part of the internal data structure to a vector shape
            result = reshape(self.dgt,  1, 4);
        end
        
        function g = eval_g(self)
            % Return the Green's function matrix g (change from Riccati parametrization to normal Green's function)
            g = ( eye(2) - self.g*self.gt ) \ ( eye(2) + self.g*self.gt );
        end
        
        function f = eval_f(self)
            % Return the anomalous Green's function matrix f (change from Riccati parametrization to normal Green's function)
            f = 2 * ( eye(2) - self.g*self.gt ) \ self.g;
        end
        
        function result = eval_ldos(self)
            result = trace(real(self.eval_g))/2;
        end
        
        function result = singlet(self)
            % Calculate the singlet component of the Green's function (proportional to iσ^y)
            f = self.eval_f;
            result = (f(1,2) - f(2,1))/2;
        end
        
        function result = triplet(self)
            % Calculate the triplet component of the Green's function (proportional to [σ^x,σ^y,σ^z]iσ^y)
            f = self.eval_f;
            result = [(f(2,2) - f(1,1))/2,                ...
                      (f(1,1) + f(2,2))/2i,               ...
                      (f(1,2) + f(2,1))/2];
        end
        
        function result = srtc(self, exchange)
            % This method returns the short-range triplet component
            % (along the exchange-field provided as an argument)
            unitvec = exchange/norm(exchange);
            result  = dot(self.triplet,unitvec) .* unitvec;
        end

        function result = lrtc(self, exchange)
            % This method returns the long-range triplet component
            % (perpendicular to the exchange-field provided as argument)
            result = self.triplet - self.srtc(exchange);
        end
    end
end