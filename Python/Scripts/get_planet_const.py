# These are Gemini written scripts. Not verified or tested yet

import numpy as np

def calculate_phi_offset(alpha_p_deg, delta_p_deg, inc_deg, raan_deg):
    # 1. Convert inputs to radians
    alpha_p = np.radians(alpha_p_deg)
    delta_p = np.radians(delta_p_deg)
    inc = np.radians(inc_deg)
    raan = np.radians(raan_deg)
    
    # 2. Compute the Planet's Spin Pole Unit Vector (N_pole)
    N_pole = np.array([
        np.cos(delta_p) * np.cos(alpha_p),
        np.cos(delta_p) * np.sin(alpha_p),
        np.sin(delta_p)
    ])
    
    # 3. Compute the Planet's Orbit Normal Unit Vector (N_orbit)
    N_orbit = np.array([
        np.sin(inc) * np.sin(raan),
        -np.sin(inc) * np.cos(raan),
        np.cos(inc)
    ])
    
    # 4. Find your GEI Frame X-Axis (Planet's True Vernal Equinox)
    V_equinox = np.cross(N_pole, N_orbit)
    X_gei = V_equinox / np.linalg.norm(V_equinox)
    
    # 5. Find the IAU J2000 Node Line X-Axis
    V_node = np.cross(np.array([0, 0, 1]), N_pole)
    X_iau = V_node / np.linalg.norm(V_node)
    
    # 6. Calculate the 1D angle between the two axes around the N_pole vector
    cos_phi = np.dot(X_iau, X_gei)
    sin_phi = np.dot(np.cross(X_iau, X_gei), N_pole)
    
    phi_rad = np.arctan2(sin_phi, cos_phi)
    return np.degrees(phi_rad) % 360.0

import numpy as np
from scipy.optimize import fsolve
import spiceypy as spice

def find_exact_equinox_jd(planet_naif_id, start_jd):
    # Load standard planetary ephemeris kernel (e.g., de440.bsp, pck00010.tpc)
    
    def subsolar_latitude_root(jd):
        # Convert Julian Date to Ephemeris Time
        et = spice.unitim(jd, 'SCLK', 'ET') 
        
        # 1. Fetch the Sun's position relative to the target planet's center 
        # in the planet's rotating body-fixed frame (e.g., 'IAU_MARS')
        state, lt = spice.spkezr('SUN', et, f'IAU_{planet_naif_id}',
                                 'NONE', planet_naif_id)
        sun_pos = state[0:3]
        
        # 2. Convert Cartesian vector to latitudinal coordinates
        r, lon, lat = spice.reclat(sun_pos)
        # This will cross zero exactly at the moment of the equinox
        return lat 
    
    # Use a numerical solver to pinpoint the exact zero crossing
    equinox_jd = fsolve(subsolar_latitude_root, start_jd)[0]
    return equinox_jd
