"""Module lambert.py contains functions for conversion between
geographic coordinates and cartesian coordinates according to
Lambert conformal conic projection - two standard parallel (2SP) case.
Implemented accroding to IOGP report, page 19: http://www.epsg.org/Portals/0/373-07-2.pdf?ver=2018-10-12-153840-577
"""
import numpy as np

class Lambert:
    """Two standard parallel (2SP) Lambert conformal conic projection
        Reference:  page 19, http://www.epsg.org/Portals/0/373-07-2.pdf?ver=2018-10-12-153840-577
    """

    # Constants
    rad_per_deg = 0.01745329251994328  # Unit coversion: Radians per degree = pi/180
    iterations  = 100                  # Number of iterations for estimating latitude in cartesian2geographic
                                       # IOGP report: "The solution should quickly converge, in 3 or 4 iterations."

    def __init__(self, semimajor_axis, inverse_flattening, standard_parallels, central_latitude, central_longitude, false_easting, false_northing):
        """Initialize ellipsoid constants and projection constants:
        Parameters:
            semimajor_axis      of the ellipsoid
            inverse_flattenging of the ellipsoid
            standard_parallels: Tuple or list containing the first and the second standard parallel
            central_latitude:   Latitude of false origin
            central_longitude:  Longitude of false origin (the Central Meridian)
            false_easting:      Easting at false origin
            false_northing:     Northing at false origin
        """
        # Ellipsoid
        self.a     = semimajor_axis                # Semi-major axis
        #b                                         # Semi-minor axis
        self.f_inv = inverse_flattening            # Inverse flattening:  f_inv = a/(a-b)
        self.f     = 1 / self.f_inv                # Flattening
        self.ecc   = np.sqrt(self.f*(2-self.f))    # Eccentricity

        # Projection constants
        self.phi_1_deg    = standard_parallels[0]  # Latitude of the first standard parallel
        self.phi_2_deg    = standard_parallels[1]  # Latitude of the second standard parallel
        self.phi_F_deg    = central_latitude       # Latitude of false origin
        self.lambda_F_deg = central_longitude      # Longitude of false origin
        self.E_F          = false_easting          # Easting at false origin
        self.N_F          = false_northing         # Northing at false origin

        # Unit conversion
        self.rad_per_deg  = Lambert.rad_per_deg
        self.phi_1_rad    = self.phi_1_deg    * self.rad_per_deg
        self.phi_2_rad    = self.phi_2_deg    * self.rad_per_deg
        self.phi_F_rad    = self.phi_F_deg    * self.rad_per_deg
        self.lambda_F_rad = self.lambda_F_deg * self.rad_per_deg


    # Macros for repeated calculations
    @staticmethod
    def _calculate_m(ecc,x): return np.cos(x) / (1 - ecc**2 * np.sin(x)**2)**0.5
    @staticmethod
    def _calculate_t(ecc,x): return np.tan(np.pi/4 - x/2) / ( (1-ecc*np.sin(x))/(1+ecc*np.sin(x)) )**(ecc/2)


    def geographic2cartesian(self, latitude, longitude):
        """Two standard parallel (2SP) Lambert conformal conic projection
        Parameters: latitude, longitude (both in degrees)
        Returns:    easting, northing (both in meters)
        Reference:  page 19, http://www.epsg.org/Portals/0/373-07-2.pdf?ver=2018-10-12-153840-577
        """
        phi_deg  = np.array(latitude)
        lmbd_deg = np.array(longitude)

        # Initialize the output vectors
        E_out = np.full(np.shape(phi_deg), np.nan)
        N_out = np.full(np.shape(phi_deg), np.nan)

        # Remove nan values from the input
        mask_ok = ~np.isnan(phi_deg) & ~np.isnan(lmbd_deg)
        phi_deg  = phi_deg[mask_ok]
        lmbd_deg = lmbd_deg[mask_ok]

        # Convert angles to radians
        lmbd_rad = lmbd_deg * self.rad_per_deg
        phi_rad  = phi_deg  * self.rad_per_deg

        # Apply all the formulas from IOGP report:
        m_1 = Lambert._calculate_m(self.ecc, self.phi_1_rad)
        m_2 = Lambert._calculate_m(self.ecc, self.phi_2_rad)

        t   = Lambert._calculate_t(self.ecc, phi_rad)
        t_1 = Lambert._calculate_t(self.ecc, self.phi_1_rad)
        t_2 = Lambert._calculate_t(self.ecc, self.phi_2_rad)
        t_F = Lambert._calculate_t(self.ecc, self.phi_F_rad)

        n = (np.log(m_1)-np.log(m_2)) / (np.log(t_1)-np.log(t_2))
        F = m_1 / (n * t_1**n)

        r   = self.a*F*(np.sign(t)*np.abs(t)**n)
        r_F = self.a*F*(np.sign(t_F)*np.abs(t_F)**n)

        theta = n*(lmbd_rad - self.lambda_F_rad)

        # Get the result
        E = self.E_F +       r*np.sin(theta)  # Easting (x-axis)
        N = self.N_F + r_F - r*np.cos(theta)  # Northing (y-axis)

        # Fill in the non-nan values to the output vectors
        E_out[mask_ok] = E
        N_out[mask_ok] = N

        return E_out, N_out


    def cartesian2geographic(self, E, N):
        """Inverse two standard parallel (2SP) Lambert conformal conic projection
        Parameters: easting, northing (both in meters)
        Returns:    latitude, longitude (both in degrees)
        Reference:  page 19, http://www.epsg.org/Portals/0/373-07-2.pdf?ver=2018-10-12-153840-577
        """
        E = np.array(E)
        N = np.array(N)

        # Initialize the output vectors
        lat_out = np.full(np.shape(E), np.nan)
        lon_out = np.full(np.shape(E), np.nan)

        # Remove nan values from the input
        mask_ok = ~np.isnan(E) & ~np.isnan(N)
        E = E[mask_ok]
        N = N[mask_ok]

        # Apply all the formulas from IOGP report
        m_1 = Lambert._calculate_m(self.ecc, self.phi_1_rad)
        m_2 = Lambert._calculate_m(self.ecc, self.phi_2_rad)

        t_1 = Lambert._calculate_t(self.ecc, self.phi_1_rad)
        t_2 = Lambert._calculate_t(self.ecc, self.phi_2_rad)
        t_F = Lambert._calculate_t(self.ecc, self.phi_F_rad)

        n = (np.log(m_1)-np.log(m_2))/(np.log(t_1)-np.log(t_2))
        F = m_1 / (n * np.sign(t_1)*np.abs(t_1)**n)
        r_F = self.a * F * (np.sign(t_F)*np.abs(t_F)**n)

        r_dash = np.sign(n) * np.sqrt(  (E-self.E_F)**2 + (r_F - (N-self.N_F))**2  )
        t_dash = np.sign(r_dash/(self.a*F)) * np.abs(r_dash/(self.a*F))**(1/n)
        theta_dash = np.arctan2( (E-self.E_F), (r_F-(N-self.N_F)) )

        # Iterative solution for the latitude
        phi = np.pi/2 - 2*np.arctan(t_dash)  # Initial value
        for _ in range(Lambert.iterations):
            phi = np.pi/2 - 2*np.arctan( t_dash * ((1 - self.ecc*np.sin(phi))/(1 + self.ecc*np.sin(phi)))**(self.ecc/2) )

        # Longitude
        lmbd = theta_dash / n + self.lambda_F_rad

        # Convert radians to degrees
        lat_out[mask_ok]  = phi  / self.rad_per_deg
        lon_out[mask_ok] = lmbd / self.rad_per_deg

        return lat_out, lon_out