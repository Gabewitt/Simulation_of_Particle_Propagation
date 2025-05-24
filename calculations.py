
import numpy as np




def random_unit_vector(forward_only=False): # Generate a random unit vector in 3D spherical space
    
    while True:
        phi = np.random.uniform(0, 2 * np.pi)
        costheta = np.random.uniform(-1, 1)
        sintheta = np.sqrt(1 - costheta**2)
        vec = np.array([sintheta * np.cos(phi), sintheta * np.sin(phi), costheta], dtype=float)
        
        if not forward_only or vec[0] > 0:
            return vec
    

def forward_biased_unit_vector(min_cos_theta=0.5): # Generate a random unit vector biased towards the +x direction in a conical shape
    
    assert 0 <= min_cos_theta <= 1.0

    while True:
        phi = np.random.uniform(0, 2 * np.pi)
        costheta = np.random.uniform(min_cos_theta, 1.0)  # +x dirction
        sintheta = np.sqrt(1 - costheta**2)
        x = costheta  
        y = sintheta * np.cos(phi)
        z = sintheta * np.sin(phi)
        vec = np.array([x, y, z])
        return vec / np.linalg.norm(vec)




def compute_absorption_probability(material, energy_ev, step_length):
   
    sigma_abs_ref = material.cross_section * material.absorption_probability
    
    E = max(energy_ev, material.reference_energy_ev)
    sigma_abs_E = sigma_abs_ref * np.sqrt(material.reference_energy_ev / E)
    return 1 - np.exp(-sigma_abs_E * step_length) # PDF of absorption over a step length



def sample_step_length(material, energy_ev=None, fast=True): # Sample distance traveled in a material before interaction
    
   
    if fast and energy_ev is not None:
        lam = material.mean_free_path(energy_ev)
    else:
        lam = material.mean_free_path(None)
    
    return np.random.exponential(lam)