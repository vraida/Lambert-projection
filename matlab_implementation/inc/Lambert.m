%Class for conversion between geographic coordinates
%and cartesian coordinates according to
%Lambert conformal conic projection - two standard parallel (2SP) case.
%Implemented accroding to IOGP report, page 19: http://www.epsg.org/Portals/0/373-07-2.pdf?ver=2018-10-12-153840-577

classdef Lambert
    %LAMBERT Two standard parallel (2SP) Lambert conformal conic projection
    %   Reference:  page 19, http://www.epsg.org/Portals/0/373-07-2.pdf?ver=2018-10-12-153840-577
    
    properties
        % Assinged in the constructor, see comments there
        a, f_inv, f, ecc, E_F, N_F
        phi_1_deg, phi_2_deg, phi_F_deg, lambda_F_deg
        phi_1_rad, phi_2_rad, phi_F_rad, lambda_F_rad
        
        % Constants
        rad_per_deg = 0.01745329251994328  % Unit coversion: Radians per degree = pi/180
        iterations  = 100                  % Number of iterations for estimating latitude in cartesian2geographic
                                           % IOGP report: "The solution should quickly converge, in 3 or 4 iterations."
    end
    
    
    methods
        function obj = Lambert(standard_parallels, central_latitude, central_longitude, false_easting, false_northing, semimajor_axis, inverse_flattening)
            %LAMBERT Construct an instance of this class
            %   Define all the projection constants.
            if  nargin > 0
                % Ellipsoid
                obj.a     = semimajor_axis;                % Semi-major axis
                %b                                         % Semi-minor axis
                obj.f_inv = inverse_flattening;            % Inverse flattening:  f_inv = a/(a-b)
                obj.f     = 1 / obj.f_inv;                 % Flattening
                obj.ecc   = sqrt(obj.f*(2-obj.f));         % Eccentricity

                % Projection constants
                obj.phi_1_deg    = standard_parallels(1);  % Latitude of the first standard parallel
                obj.phi_2_deg    = standard_parallels(2);  % Latitude of the second standard parallel
                obj.phi_F_deg    = central_latitude;       % Latitude of false origin
                obj.lambda_F_deg = central_longitude;      % Longitude of false origin
                obj.E_F          = false_easting;          % Easting at false origin
                obj.N_F          = false_northing;         % Northing at false origin

                % Unit conversion
                obj.phi_1_rad    = obj.phi_1_deg    * obj.rad_per_deg;
                obj.phi_2_rad    = obj.phi_2_deg    * obj.rad_per_deg;
                obj.phi_F_rad    = obj.phi_F_deg    * obj.rad_per_deg;
                obj.lambda_F_rad = obj.lambda_F_deg * obj.rad_per_deg;
            end
        end
        
        
        function [E_out, N_out] = geographic2cartesian(obj, latitude, longitude)
            %GEOGRAPHIC2CARTESIAN Two standard parallel (2SP) Lambert conformal conic projection
            %   Parameters: latitude, longitude (both in degrees)
            %   Returns:    easting, northing (both in meters)
            %   Reference:  page 19, http://www.epsg.org/Portals/0/373-07-2.pdf?ver=2018-10-12-153840-577
            phi_deg  = latitude;
            lmbd_deg = longitude;

            % Initialize the output vectors
            E_out = nan(size(phi_deg));
            N_out = nan(size(phi_deg));

            % Remove nan values from the input
            mask_ok = ~isnan(phi_deg) & ~isnan(lmbd_deg);
            phi_deg  = phi_deg(mask_ok);
            lmbd_deg = lmbd_deg(mask_ok);

            % Convert angles to radians
            lmbd_rad = lmbd_deg * obj.rad_per_deg;
            phi_rad  = phi_deg  * obj.rad_per_deg;

            % Apply all the formulas from IOGP report:
            m_1 = Lambert.calculate_m(obj.ecc, obj.phi_1_rad);
            m_2 = Lambert.calculate_m(obj.ecc, obj.phi_2_rad);

            t   = Lambert.calculate_t(obj.ecc, phi_rad);
            t_1 = Lambert.calculate_t(obj.ecc, obj.phi_1_rad);
            t_2 = Lambert.calculate_t(obj.ecc, obj.phi_2_rad);
            t_F = Lambert.calculate_t(obj.ecc, obj.phi_F_rad);

            n = (log(m_1)-log(m_2)) / (log(t_1)-log(t_2));
            F = m_1 / (n * t_1^n);

            r   = obj.a*F*(t.^n);
            r_F = obj.a*F*(t_F^n);

            theta = n*(lmbd_rad - obj.lambda_F_rad);

            % Get the result
            E = obj.E_F +       r.*sin(theta);  % Easting (x-axis)
            N = obj.N_F + r_F - r.*cos(theta);  % Northing (y-axis)

            % Fill in the non-nan values to the output vectors
            E_out(mask_ok) = E;
            N_out(mask_ok) = N;
        end
        
        
        function [lat_out, lon_out] = cartesian2geographic(obj, E, N)
            %CARTESIAN2GEOGRAHPIC Inverse two standard parallel (2SP) Lambert conformal conic projection
            %Parameters: easting, northing (both in meters)
            %Returns:    latitude, longitude (both in degrees)
            %Reference:  page 19, http://www.epsg.org/Portals/0/373-07-2.pdf?ver=2018-10-12-153840-577

            % Initialize the output vectors
            lat_out = nan(size(E));
            lon_out = nan(size(E));

            % Remove nan values from the input
            mask_ok = ~isnan(E) & ~isnan(N);
            E = E(mask_ok);
            N = N(mask_ok);

            % Apply all the formulas from IOGP report
            m_1 = Lambert.calculate_m(obj.ecc, obj.phi_1_rad);
            m_2 = Lambert.calculate_m(obj.ecc, obj.phi_2_rad);

            t_1 = Lambert.calculate_t(obj.ecc, obj.phi_1_rad);
            t_2 = Lambert.calculate_t(obj.ecc, obj.phi_2_rad);
            t_F = Lambert.calculate_t(obj.ecc, obj.phi_F_rad);

            n = (log(m_1)-log(m_2))/(log(t_1)-log(t_2));
            F = m_1 / (n * t_1^n);
            r_F = obj.a * F * (t_F^n);

            r_dash = sign(n) * sqrt(  (E-obj.E_F).^2 + (r_F - (N-obj.N_F)).^2  );
            t_dash = sign(r_dash/(obj.a*F)) .* abs(r_dash/(obj.a*F)).^(1/n);
            theta_dash = atan2( (E-obj.E_F), (r_F-(N-obj.N_F)) );

            % Iterative solution for the latitude
            phi = pi/2 - 2*atan(t_dash);  % Initial value
            for ii = 1 : obj.iterations
                phi = pi/2 - 2*atan( t_dash .* ((1 - obj.ecc*sin(phi))./(1 + obj.ecc*sin(phi))).^(obj.ecc/2) );
            end
            
            % Longitude
            lmbd = theta_dash / n + obj.lambda_F_rad;

            % Convert radians to degrees
            lat_out(mask_ok) = phi  / obj.rad_per_deg;
            lon_out(mask_ok) = lmbd / obj.rad_per_deg;
        end
        
    end
    
    
    methods(Static)
        % Macros for repeated calculations
        function m = calculate_m(ecc,x)
            m = cos(x) / (1 - ecc^2 * sin(x)^2)^0.5;
        end
        function t = calculate_t(ecc,x)
            t = tan(pi/4 - x/2) ./ ( (1-ecc*sin(x))./(1+ecc*sin(x)) ).^(ecc/2);
        end
    end
    
end

