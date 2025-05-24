import numpy as np

class Material:
    def __init__(self, name, cross_section, density, length_x, absorption_probability, mass_number):
        self.name = name
        self.cross_section = cross_section  
        self.density = density              
        self.length_x = length_x            
        self.absorption_probability = absorption_probability
        self.mass_number = mass_number      
        
        A = self.mass_number
        self.energy_loss_factor = (A**2 + 1) / (A + 1)**2

        self.reference_energy_ev = 0.025  # Reference energy in eV (thermal neutron energy at room temperature)

    def mean_free_path(self, energy_ev=None):
        
        if energy_ev is None:
            return 1 / self.cross_section

        scaling = np.sqrt(energy_ev / self.reference_energy_ev)
        return scaling / self.cross_section #energy dependent mean free path


class Particle:
    def __init__(self, name, mass, charge, initial_velocity):
    
        self.name = name
        self.mass = mass  
        self.charge = charge  
        self.initial_velocity = initial_velocity 
        self.energy = 0.5 * self.mass * np.linalg.norm(initial_velocity)**2 
        self.energy_ev = self.energy / 1.602e-19  
    
    def is_thermal(self, temperature):
        kb = 8.617e-5  
        Eth = 3/2 * kb * temperature
        
        return self.energy_ev <= Eth

    def set_thermal_velocity(self, temperature):
        kb_J = 1.380649e-23  
        E_thermal = 1.5 * kb_J * temperature 
        v = (2 * E_thermal / self.mass)**0.5  

        self.initial_velocity = [v, 0, 0]
        self.energy = E_thermal
        self.energy_ev = E_thermal / 1.602e-19
        
        
        
        


        