classdef quaternion2 < rotation
    %QUATERNION CLASS: 
    %
    % This class represents the quaternion2 of rotation following the
    % convention in Shuster 1993:
    %
    %  - The scalar part is in the last quaternion2 component;
    %  - Multiplication is implemented with the left convention 
    % (see Shuster 1993,eq 172). 
    %
    % The order of  quaternion2 multiplication is the same of the rotation 
    % matrix composition.
    %
    %
    % Author: Gabriele Ligorio
    % Version: 1.0
    % Date: 04-11-2015
    %
    
    methods
        
        %% Constructor
        function obj = quaternion2(varargin)
            
            if nargin == 0
                
                obj.value = [0; 0; 0; 1];
                obj.len = 1;
                
            elseif nargin == 1
                % Construct from another object (Copy Constructor)
                if isa(varargin{1},'quaternion2')
                    obj.value = varargin{1}.getValue();
                    obj.len = varargin{1}.getN();
                    
                elseif isa(varargin{1},'rotMatrix')
                    R = varargin{1}.getValue();
                    obj = quaternion2(RM2quat(R));
                    
                elseif isa(varargin{1},'EA')
                    ea = varargin{1}.getValue();
                    eaSeq = varargin{1}.getSeq();
                    obj = quaternion2(RM2quat(ea2RM(ea,eaSeq)));
                    
                elseif length(varargin{1}) == 1
                     % Construct quaternion2 by specifying only the length
                    obj.value = [zeros(3,varargin{1}); ones(1,varargin{1})];
                    obj.len = varargin{1};
                    
                else
                    % Construct quaternion2 by passing a 4xN matrix
                    if (size(varargin{1},1) ~= 4)
                        error('Input Quaternion must be 4xN')
                    else
                        obj.value = varargin{1};
                        obj.len = size(obj.value,2);
                    end
                end
            elseif nargin == 2
                    % Construct quaternion2 from a Nx4 matrix
                    if (size(varargin{1},2) == 4) && (strcmp(varargin{2},'conv'))
                        q = varargin{1};
                        obj.value = [q(:,2:4) q(:,1)]';
                        obj.len = size(obj.value,2);
                    else
                        error('Check input parameters')
                    end
            else
                error('Too many input arguments');
            end
     
        end
        
        %% Get Methods 
        
        % Get quaternion2 value
        function rslt = getValue(obj,pos)
            
            if nargin == 1
                rslt = obj.value;
            elseif nargin == 2
                rslt = obj.value(:,pos);
            else
                error('Too many input parameters');
            end            
            
        end
        
        % Get quaternion2 dimension
        function rslt = getN(obj)
            rslt = obj.len;
        end
        
        % Get Scalar part
        function rslt = getScalar(obj)
            rslt = obj.value(4,:);
        end
        
        %% Set methods
        function setValue(obj,in,pos)
            
            if size(in,1)== 4 || isa(in,'quaternion2')
                
                if nargin == 2
                    obj.value = in;
                    obj.len = size(in,2);
                elseif nargin == 3
                    obj.value(:,pos) = in;
                else
                    error('Too many input parameters');
                end
                
            else
                error('Wrong input size or type')
            end
            
        end
        
        %% Inherited methods
        
        % Select quaternion2 values at x
        function rslt = select(obj,x)
            %temp = obj.getValue;
            rslt = quaternion2(obj.getValue(x));
        end
        
        % Return the conjugate quaternion2 object
        function rslt = inv(obj)
            q = obj.getValue;
            rslt = quaternion2([-q(1:3,:); q(4,:)]);            
        end
                              
        % Brute-force quaternion2 normalization
        function unit(obj)
            obj.setValue(obj.getValue./repmat(obj.norm(),4,1));
        end
        
        % Replicate quaternion2 value
        function rslt = rep(obj,N)
            q = obj.getValue();
            q = repmat(q,1,N);
            rslt = quaternion2(q);
        end
        
        % Vector rotation
        function rslt = rotate(obj,vect)
            if (size(vect,1)~= 3)
                error('Wrong input dimension!')
            else
                q1 = obj.getValue();
                q2 = vect;
                % Calculate vector portion of quaternion2 product
                % vec = s2*v1 + s1*v2 + cross(v2,v1) with s2 = 0!!
                vec12 = [q1(4,:).*q2(1,:); q1(4,:).*q2(2,:); q1(4,:).*q2(3,:)]+...
                    [ q2(2,:).*q1(3,:)-q2(3,:).*q1(2,:); ...
                    q2(3,:).*q1(1,:)-q2(1,:).*q1(3,:); ...
                    q2(1,:).*q1(2,:)-q2(2,:).*q1(1,:)];
                
                % Calculate scalar portion of quaternion2 product
                % scalar = s1*s2 - dot(v1,v2) with s2 = 0!
                scalar12 = - q2(1,:).*q1(1,:) - q2(2,:).*q1(2,:) - q2(3,:).*q1(3,:);
                
                % q3 = q*[vect;0]
                q3 = [vec12; scalar12];
                
                % Now multiply calculate only the vectorial part of q3*q1.inv()
                % q1 signs are inverted instead of calling q1.inv()
                rslt = [q1(4,:).*q3(1,:); q1(4,:).*q3(2,:); q1(4,:).*q3(3,:)] + ...
                    [q3(4,:).*-q1(1,:); q3(4,:).*-q1(2,:); q3(4,:).*-q1(3,:)]+...
                    [ -q1(2,:).*q3(3,:)+q1(3,:).*q3(2,:); ...
                    -q1(3,:).*q3(1,:)+q1(1,:).*q3(3,:); ...
                    -q1(1,:).*q3(2,:)+q1(2,:).*q3(1,:)];
                
            end
            
        end
        
        % Mean rotation in an time window
        function rslt = average(obj,win)
            % Time window of interest
            temp = obj.getValue(win(1):win(2));
            
            % Mean quaternion2
            rslt = quaternion2(mean(temp,2));
            
            % Enforce unitary norm constraint
            rslt.unit;
        
        end
        
        function proj(obj,w_b,T,init)
            N = length(w_b);
            q = zeros(4,N);
            
            if nargin == 4
                if isa(init,'quaternion2')
                    q(:,1) = init.getValue();
                else
                    error('Initial condition: wrong type!');
                end
            else
                q(:,1) = [0 0 0 1]';
            end
            
            for i = 1:N-1                
                
                Phi = quatProjMat(w_b(:,i),T);
                q(:,i+1) = Phi*q(:,i);
                
            end
            
            obj.setValue(q);
        end
        
        %% Quaternion methods
        
        % Get quaternion2 value with matlab convention
        function rslt = getConvertedValue(obj,pos)
            
            if nargin == 1
                rslt = obj.value;
            elseif nargin == 2
                rslt = obj.value(:,pos);
            else
                error('Too many input parameters');
            end           
            
            rslt = [rslt(4,:); rslt(1:3,:)]';
        end
        
        
        % Return the quaternion2 norm (it should be 1)
        function rslt = norm(obj)
            rslt = sqrt(sum(obj.getValue.^2));
        end
        
        % Xi operator (Shuster 1993, eq 307)
        function rslt = xi(obj)            
            rslt = zeros(4,3,obj.getN);
            
            for i = 1:obj.getN
            q = obj.getValue(i);
            
            rslt(:,:,i) = [ q(4)  -q(3)  q(2);
                            q(3)   q(4) -q(1);
                           -q(2)   q(1)  q(4);
                           -q(1)  -q(2) -q(3)];
            end
        end
        
        % Multipication operator (Shuster 1993, eq 308)
        function rslt = rightMultOper(obj)
            q = obj.getValue();
            q = reshape(q,4,1,obj.getN);
            rslt = [obj.xi q];
        end
        
                
        % Psi operator (Shuster 1993, eq 317)
        function rslt = psi(obj)
                        
            rslt = zeros(4,3,obj.getN);
            
            for i = 1:obj.getN
                q = obj.getValue(i);
            
                rslt(:,:,i) = [  q(4)   q(3)  -q(2);
                                -q(3)   q(4)   q(1);
                                 q(2)  -q(1)   q(4);
                                -q(1)  -q(2)  -q(3)];
            end            

        end
        
        % Multipication operator (Shustet 1993, eq 316)
        function rslt = leftMultOper(obj)
            q = obj.getValue();
            q = reshape(q,4,1,obj.getN);
            rslt = [obj.psi q];
        end
        
        %% Operators Overloading 
               
        % Quaternion sum
        function objSum = plus(obj1,obj2)
            
                q2 = obj2.getValue;
                q1 = obj1.getValue;
                
                objSum = quaternion2(q1+q2);
                
        end
        
        % Quaternion difference
        function objMin = minus(obj1,obj2)
            
            q1 = obj1.getValue;
            q2 = obj2.getValue;
            
            objMin = quaternion2(q1-q2);
            
        end
        
        % Quaternion multiplication, left convention (see Shuster 1993, eq 172)
        function objProd = mtimes(obj1,obj2)
            
            q1 = obj1.getValue;
            q2 = obj2.getValue;
            
            % Calculate vector portion of quaternion2 product
            % vec = s2*v1 + s1*v2 + cross(v2,v1)
            vec = [q2(4,:).*q1(1,:); q2(4,:).*q1(2,:); q2(4,:).*q1(3,:)] + ...
                [q1(4,:).*q2(1,:); q1(4,:).*q2(2,:); q1(4,:).*q2(3,:)]+...
                [ q2(2,:).*q1(3,:)-q2(3,:).*q1(2,:); ...
                q2(3,:).*q1(1,:)-q2(1,:).*q1(3,:); ...
                q2(1,:).*q1(2,:)-q2(2,:).*q1(1,:)];
            
            % Calculate scalar portion of quaternion2 product
            % scalar = s1*s2 - dot(v1,v2)
            scalar = q2(4,:).*q1(4,:) - q2(1,:).*q1(1,:) - ...
                q2(2,:).*q1(2,:) - q2(3,:).*q1(3,:);            
            
            objProd = quaternion2([vec; scalar]);
        end
        
    end
    
end

