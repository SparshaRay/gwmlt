"""
Decihertz Detector Classes
"""

import numpy as np
from numpy import cos, sin, pi

from astropy.time import Time
from astropy import coordinates as coord
from astropy import units as u

from pycbc.types.timeseries import TimeSeries

from gwmlt.constants import C_SI, AU_SI


class IndIGO_D :

    """
    Detector class for the Heliocentric configuration of the IndIGO-D detector.
    Duck-typed to be compatible with the `project_wave` method of PyCBC `Detector` class.
    Assumes stationary detector and long wavelength approximation (LWA) for the projection.
    Based on : https://arxiv.org/pdf/2601.06956
    """

    def __init__(self, 
        name: str = "IndIGO-D", 
        lead_offset: float = 20.0, 
        L_SI: float = 1e6
    ) -> None :

        """
        Create a new instance of the IndIGO-D detector.

        Parameters
        ----------
        name : str, optional
            The name of the detector. Default is "IndIGO-D".
        lead_offset : float, optional
            Angle by which the constellation leads earth in degrees. Default is 20.0.
        L_SI : float, optional
            The arm length of the detector in meters. Default is 1e6 (1e3 km).
        """

        self.name = name
        self.lead_offset = lead_offset
        self.L_SI = L_SI


    def _get_sc_loc(self, gps_time) :

        """
        Get the position of the 3 spacecrafts in heliocentric frame.
        """

        a, L  = AU_SI, self.L_SI
        e_ecc = L / (2.0 * a)
        alpha = e_ecc
        epsilon = np.sqrt(3) * alpha

        earth_bary = coord.get_body_barycentric('earth', Time(gps_time, format='gps'))
        earth_ecliptic = coord.SkyCoord(earth_bary, frame='icrs').barycentricmeanecliptic

        earth_phi = earth_ecliptic.lon.rad 
        phi_sc = earth_phi + np.radians(self.lead_offset)

        psi2 = (phi_sc       ) - e_ecc * sin(phi_sc       )
        psi3 = (phi_sc - pi/2) - e_ecc * sin(phi_sc - pi/2)
        
        # SC1 (Vertex)
        x1 = a * cos(phi_sc)
        y1 = a * sin(phi_sc)
        z1 = 0.0
        
        # SC2
        x2 = + a * (cos(psi2) + e_ecc) * cos(epsilon)
        y2 = + a * np.sqrt(1 - e_ecc**2) * sin(psi2)
        z2 = + a * (cos(psi2) + e_ecc) * sin(epsilon)

        # SC3
        x3 = - a * np.sqrt(1 - e_ecc**2) * sin(psi3)
        y3 = + a * (cos(psi3) + e_ecc) * cos(epsilon)
        z3 = + a * (cos(psi3) + e_ecc) * sin(epsilon)

        return np.array([x1, y1, z1]), np.array([x2, y2, z2]), np.array([x3, y3, z3])
    

    def _reference_pol_tensors(self, longitude, latitude) :

        """
        Get the reference polarization tensors for a given source position in ecliptic coordinates.
        """

        u_vec = np.array([sin(longitude), -cos(longitude), 0])
        v_vec = np.array([-sin(latitude)*cos(longitude), -sin(latitude)*sin(longitude), cos(latitude)])
        
        e_plus  = np.outer(u_vec, u_vec) - np.outer(v_vec, v_vec)
        e_cross = np.outer(u_vec, v_vec) + np.outer(v_vec, u_vec)
        return e_plus, e_cross


    def _antenna_pattern(self, ra, dec, polarization, gps_time) :

        """
        Get the Fp Fc for a given source position and polarization angle at a given time.
        """

        # Convert ICRS RA/Dec to Barycentric True Ecliptic
        source_c = coord.SkyCoord(ra=ra*u.rad, dec=dec*u.rad, frame='icrs')
        ecliptic = source_c.barycentricmeanecliptic
        lon, lat = ecliptic.lon.rad, ecliptic.lat.rad
        
        sc1, sc2, sc3 = self._get_sc_loc(gps_time)

        # Arm vectors
        link12 = (sc2 - sc1) / np.linalg.norm(sc2 - sc1)
        link13 = (sc3 - sc1) / np.linalg.norm(sc3 - sc1)

        d12 = np.outer(link12, link12)
        d13 = np.outer(link13, link13)

        e_plus, e_cross = self._reference_pol_tensors(lon, lat)

        # Rotate to detector's polarization frame
        E_plus  = + cos(2*polarization)*e_plus + sin(2*polarization)*e_cross
        E_cross = - sin(2*polarization)*e_plus + cos(2*polarization)*e_cross
        
        Fp = 0.5 * np.trace(np.dot(d12 - d13, E_plus))
        Fc = 0.5 * np.trace(np.dot(d12 - d13, E_cross))
        
        return Fp, Fc
    

    def time_delay_from_earth_center(
        self, 
        ra:float, 
        dec:float, 
        gps_time:float
    ) -> float :

        """
        Time delay (in seconds) between geocenter and IndIGO-D vertex.
        The maximum possible time delay is ~ 173 seconds.

        Parameters
        ----------
        ra : float
            Right Ascension of the source in radians.
        dec : float
            Declination of the source in radians.
        gps_time : float
            GPS time at which to compute the time delay.

        Returns
        -------
        float
            Time delay in seconds between geocenter and IndIGO-D vertex.
        """

        source_c = coord.SkyCoord(ra=ra*u.rad, dec=dec*u.rad, frame='icrs')
        ecliptic = source_c.barycentricmeanecliptic
        lon, lat = ecliptic.lon.rad, ecliptic.lat.rad
        
        # GW propagation direction
        k_hat = np.array([-cos(lon)*cos(lat), -sin(lon)*cos(lat), -sin(lat)])
        
        # Earth's position (barycentric frame)
        earth_bary = coord.get_body_barycentric('earth', Time(gps_time, format='gps'))
        earth_ecliptic = coord.SkyCoord(earth_bary, frame='icrs').barycentricmeanecliptic
        r_earth = np.array([earth_ecliptic.cartesian.x.to('m').value,
                            earth_ecliptic.cartesian.y.to('m').value,
                            earth_ecliptic.cartesian.z.to('m').value])

        # Detector position
        r_det, _, _ = self._get_sc_loc(gps_time)
        
        # Project baseline onto k_hat
        dt = np.dot(r_det - r_earth, k_hat) / C_SI
        return dt


    def project_wave(
        self, 
        hp: TimeSeries, 
        hc: TimeSeries, 
        ra: float, 
        dec: float, 
        polarization: float
    ) -> TimeSeries :

        """
        Project hp hc onto the detector to get strain.
        This method is duck-typed to be compatible with the `project_wave` method of PyCBC `Detector` class.

        Parameters
        ----------
        hp : TimeSeries
            The h_plus polarization timeseries.
        hc : TimeSeries
            The h_cross polarization timeseries.
        ra : float
            Right Ascension of the source in radians.
        dec : float
            Declination of the source in radians.
        polarization : float
            Polarization angle of the source in radians.

        Returns
        -------
        TimeSeries
            The projected strain timeseries for the detector.
        """

        # Reference time is at the end of the waveform
        t_eval = float(hp.end_time)
        
        # Get antenna patterns
        Fp, Fc = self._antenna_pattern(ra, dec, polarization, t_eval)
        
        # Project the time series
        strain = hp * Fp + hc * Fc
        
        # Offset by geocenter-detector time delay
        dt = self.time_delay_from_earth_center(ra, dec, t_eval)
        strain.start_time += dt
        
        return strain