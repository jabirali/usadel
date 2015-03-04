% Define a data structure to describe the state of the physical system
% for a given position and energy. This is done by describing the
% Riccati parametrized Green's function 'g', it's tilde conjugate 'gt', 
% and their first derivatives 'dg' and 'dgt' for that configuration.
%
% This class is mainly intended for use with differential equation
% solvers, and therefore provides the method 'vectorize' to pack the
% internal variables in a vector format, and constructor State(...)
% that is able to unpack this vector format. Alternatively, the
% vectorization can be performed without instantiating the class at all,
% by calling the static methods 'pack' and 'unpack'.
%
% Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
% Created 2015-02-14
% Updated 2015-02-24


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
            %  (i)   four 2x2 matrices, which should correspond to the
            %        matrices g, dg, gt, dgt, in that order;
            %  (ii)  one 16-element complex vector, which should be produced
            %        by either the 'vectorize' metohd or the 'pack' method;
            %  (iii) no arguments, i.e. the empty constructor.
            switch nargin
                case 1
                    % If we get one input, then assume that we got a vector
                    % created by 'vectorize', and reverse the procedure
                    
                    [self.g,self.dg,self.gt,self.dgt] = self.unpack(varargin{1});
                    
                case 4
                    % If we get four input arguments, then assume that
                    % these correspond to g, dg, gt, and dgt, respectively.
                    % Multiply with a 2x2 identity matrix in case the
                    % input arguments were scalars and not matrices.
                    
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
            result = self.pack(self.g,self.dg,self.gt,self.dgt);
        end
                
        function g = eval_g(self)
            % Return the Green's function matrix g, i.e. convert the
            % Riccati parameter self.g to a normal Green's function.
            g = ( eye(2) - self.g*self.gt ) \ ( eye(2) + self.g*self.gt );
        end
        
        function f = eval_f(self)
            % Return the anomalous Green's function matrix f, i.e. convert
            % the Riccati parameter self.g to a normal Green's function.
            f = ( eye(2) - self.g*self.gt )  \ (2 * self.g);
        end
        
        function result = eval_ldos(self)
            result = trace(real(self.eval_g))/2;
        end
        
        function result = singlet(self)
            % Calculate the singlet component of the Green's function,
            % i.e. the component proportional to iσ^y.
            f = self.eval_f;
            result = (f(1,2) - f(2,1))/2;
        end
        
        function result = triplet(self)
            % Calculate the triplet component of the Green's function,
            % i.e. the component proportional to [σ^x,σ^y,σ^z] iσ^y.
            f = self.eval_f;
            result = [(f(2,2) - f(1,1))/2,                ...
                      (f(1,1) + f(2,2))/2i,               ...
                      (f(1,2) + f(2,1))/2];
        end
        
        function result = srtc(self, exchange)
            % This method returns the short-range triplet component, i.e.
            % the triplet component *along* the exchange field, where the
            % exchange field should be provided as an argument.
            unitvec = exchange/norm(exchange);
            result  = dot(unitvec,self.triplet) .* unitvec;
        end

        function result = lrtc(self, exchange)
            % This method returns the long-range triplet component, i.e.
            % the triplet component *perpendicular* to the exchange field,
            % where the exchange field should be provided as argument.
            result = self.triplet - self.srtc(exchange);
        end
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define static methods (available without object instantiation)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    methods (Static)
        function [g,dg,gt,dgt] = unpack(vector)
            % This method is used to unpack a 16 element state vector into
            % the Riccati parametrized Green's function that it represents.
            args = reshape(vector, 2, 8);
            g    = args(:,1:2);
            dg   = args(:,3:4);
            gt   = args(:,5:6);
            dgt  = args(:,7:8);
        end
        
        function vector = pack(g,dg,gt,dgt)
            % This method is used to pack Riccati parametrized Green's
            % functions into a 16 element complex state vector.
            vector = reshape([g dg gt dgt], 16, 1);
        end
    end
end