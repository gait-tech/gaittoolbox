classdef dipole
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % {n} wrt {l}
        n_l;
        qln;
        % {l} wrt {n}
        l_n;
        %% Loop parameters
        M = 0.15;                   % z-oriented, A*m^2
        %% Air
        mu = 4*pi*1e-7;             % H/m (oppure N/(A^2))
    end
    
    methods
        %% Constructor
        function obj = dipole(qln, l_n, M)
            if isa(qln,'quaternion2')
                obj.qln = qln;
            else
                error('dipole@dipole: Wrong input 1.');
            end
                                  
            if ~isequal(size(l_n), [3 1])
                error('dipole@dipole: Wrong input 2.');
            end
            % magnetic moment, M
            if nargin > 2
                if isscalar(M)
                    obj.M = M;
                else
                    error('dipole@dipole: Wrong input 3.');
                end
            end
            % {n} position in {l} - position of the dipoles origin in the
            % local co-oridinate frame
            obj.n_l = -obj.qln.rotate(l_n);
            % {l} position in {n} - position of the dipole's origin in the
            % navigation frame
            obj.l_n = l_n;
        end
        
        %% Magnetic field calculation
        % p_n is the point (in the navigation frame {n}) where I want to
        % evaluate the magnetic field
        function rslt = calc_B(obj,p_n)
            % point in the navigation frame expressed relative to the
            % dipole's origin
            p_l = obj.n_l + obj.qln.rotate(p_n);
            % Cylindric coordinates
            % the x/y component of the dipole with respect to the cartesian
            % coordinates defining the magnetic moment, see Fig. 5, this
            % will vary depending on the quaternion
            phi = atan2(p_l(1),p_l(2)); 
            % horizontal distance between the observation point and the
            % origin of the dipole
            r = norm(p_l(1:2)); 
            % vertical distance between the observation point and the
            % origin of the dipole
            z = p_l(3);
            
            % N/(A*m^3) = mu_a * M / (4 * pi * (r^2 + z^2)^(5/2))
            cost_vector =  obj.mu*obj.M/(4*pi*(r^2+z^2)^(5/2)); 
            
            % microTesla - Equation (5)
            Bx = 3*cost_vector*r*z*sin(phi)*1e6;                 
            By = 3*cost_vector*r*z*cos(phi)*1e6;                 
            Bz = cost_vector.*(2*z.^2 - r.^2)*1e6;               
            
            rslt = [Bx;By;Bz];
             
        end
           
        % Calculate and show magnetic field
        function show_B(obj)
            rx = 0.1; ry = 0.1; rz = 0.1;
%             rx = 0.05; ry = 0.05; rz = 0.05;
            res = 0.01;
            
            points_l = [-rx:res:rx; -ry:res:ry; -rz:res:rz];

            points_n = obj.qln.inv().rotate(points_l - repmat(obj.n_l,1,length(points_l)));

            x = points_n(1,:);
            y = points_n(2,:);
            z = points_n(3,:);
            
            Nx = length(x);
            Ny = length(y);
            Nz = length(z);

            [X, Y, Z] = meshgrid(x,y,z);
            
            Bx = zeros(size(X));
            By = zeros(size(X));
            Bz = zeros(size(X));
            
            steps = Nx*Ny*Nz;
            
            h = waitbar(0,'Magnetic field calculation...');
            for i = 1:Nx
                for j = 1:Ny
                    for k = 1:Nz
                        p = [X(i,j,k) Y(i,j,k) Z(i,j,k)]';
                        rslt = obj.calc_B(p);
                        Bx(i,j,k) = rslt(1);
                        By(i,j,k) = rslt(2);
                        Bz(i,j,k) = rslt(3);
                        waitbar((i*Nx*Ny + j*Ny + k)/steps);
                    end
                end
            end
            
            Bx(isnan(Bx)) = 0;
            By(isnan(By)) = 0;
            Bz(isnan(Bz)) = 0;
            
            close(h);
            bb = 0.06;
            
            
            plot_field = streamslice(X,Y,Z,Bx,By,Bz,obj.l_n(1),obj.l_n(2),obj.l_n(3));
            set(plot_field,'Color','red');
            hold on
            quiver3(X,Y,Z,Bx,By,Bz,'linewidth',2);
            grid on;
            xlabel('x [m]'), ylabel('y [m]'), zlabel('z [m]')
            view(gca,[-193 18])
            line([(obj.l_n(1)-bb) (obj.l_n(1)+bb)],[obj.l_n(2) obj.l_n(2)],[obj.l_n(3) obj.l_n(3)],'color','black')
            line([obj.l_n(1) obj.l_n(1)],[(obj.l_n(2)-bb) (obj.l_n(2)+bb)],[obj.l_n(3) obj.l_n(3)],'color','black')
            line([obj.l_n(1) obj.l_n(1)],[obj.l_n(2) obj.l_n(2)],[(obj.l_n(3)-bb) (obj.l_n(3)+bb)],'color','black')
        end
        
    end
    
end

