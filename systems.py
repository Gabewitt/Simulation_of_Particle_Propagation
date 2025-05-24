import numpy as np
from constants import Material, Particle
from properties import (
    Water, Graphite, Lead, Air, LightWater)
from calculations import (
random_unit_vector,
forward_biased_unit_vector,
compute_absorption_probability,
sample_step_length
)




def propagate_particle(
    material, 
    particle, 
    x_entry=0, 
    x_start=-50,
    max_steps=1000,
    thermal_mode=True,
    One_Directional_Propagation=False,
    min_cos_theta=0.5,
    Fast_Ballistic_Propagation=True):

    # Initialize position and path
    pos = np.array([x_start, 0, 0])
    path = [pos.copy().tolist()]
    energies = [particle.energy_ev]
    energies_in_material = []
    end_position = None

    # Straight line entry steps
    n_entry_steps = 10
    for i in range(1, n_entry_steps + 1):
        xi = x_start + (x_entry - x_start) * (i / n_entry_steps)
        pos = np.array([xi, 0.0, 0.0])
        path.append(pos.copy().tolist())
        energies.append(particle.energy_ev)

    for _ in range(max_steps):
        
        if thermal_mode: #if true, particle is thermal 
            direction = random_unit_vector(forward_only=One_Directional_Propagation)
        else:
            direction = forward_biased_unit_vector(min_cos_theta)

        step = sample_step_length(
            material,
            energy_ev=particle.energy_ev if not thermal_mode else None,
            fast=Fast_Ballistic_Propagation
        )

        pos = pos + step * direction
        path.append(pos.copy().tolist())
        energies.append(particle.energy_ev)

        if pos[0] > 0:
            energies_in_material.append(particle.energy_ev)

        
        if not thermal_mode: # if particle is not thermal
            
            particle.energy_ev *= material.energy_loss_factor
            particle.energy   = particle.energy_ev * 1.602e-19  

        p_abs = compute_absorption_probability(material,
                                       particle.energy_ev,
                                       step) 
        # Check for absorption
        if np.random.rand() < p_abs:
            print(f"Absorbed particle at x = {pos[0]:.3f} m, energy = {particle.energy_ev:.2f} eV")
            end_position = pos[0]
            break
        
        # Check if particle exited the slab
        if not (0 <= pos[0] <= material.length_x):
            end_position = pos[0]
            print(f"Neutron exited slab at x = {pos[0]:.3f} m (y = {pos[1]:.3f}, z = {pos[2]:.3f})")
            break

    path_with_energy = [p + [e] for p, e in zip(path, energies)]

    return path, path_with_energy, None, energies, energies_in_material, end_position






def propagate_moderated_particle(
    material,
    particle,
    x_entry=0,
    x_start=-50,
    max_steps=1000,
    temperature=293,
    One_Directional_Propagation=False,
    min_cos_theta=0.5,
    Fast_Ballistic_Propagation=True
):
    
    pos = np.array([x_start, 0, 0])
    path = [pos.copy().tolist()]
    energies = [particle.energy_ev]
    thermalized_index = None
    energies_in_material = []
    end_position = None

    
    for i in range(1, 11):
        xi = x_start + (x_entry - x_start) * (i / 10)
        pos = np.array([xi, 0.0, 0.0])
        path.append(pos.tolist())
        energies.append(particle.energy_ev)


    for _ in range(max_steps):
        
        if particle.is_thermal(temperature): 
            
            
            if thermalized_index is None:
                thermalized_index = len(path) - 1

            
            direction = random_unit_vector(forward_only=One_Directional_Propagation)
            
            
            step = sample_step_length(material, energy_ev=None, fast=False)

        else:
            
            
            direction = forward_biased_unit_vector(min_cos_theta)
            step = sample_step_length(
                material,
                energy_ev=particle.energy_ev,
                fast=Fast_Ballistic_Propagation
            )

        
        pos = pos + step * direction
        path.append(pos.tolist())
        energies.append(particle.energy_ev)

        if pos[0] > 0:
            energies_in_material.append(particle.energy_ev)

        
        p_abs = compute_absorption_probability(material, particle.energy_ev, step)
        
        # Check for absorption
        if np.random.rand() < p_abs:
            print(f"Absorbed particle at x = {pos[0]:.3f} m, energy = {particle.energy_ev:.2f} eV")
            end_position = pos[0]
            break

        
        if not particle.is_thermal(temperature):
             
            particle.energy_ev *= material.energy_loss_factor
            particle.energy   = particle.energy_ev * 1.602e-19  

            
        # Check if particle exited the slab
        if not (0 <= pos[0] <= material.length_x):
            print(f"Neutron exited slab at x = {pos[0]:.3f} m")
            end_position = pos[0]
            break

    
    path_with_energy = [p + [e] for p, e in zip(path, energies)]
    return path, path_with_energy, thermalized_index, energies, energies_in_material, end_position



def simulate_multiple_neutrons(
    material,
    particle,
    simulation_mode,           
    n_particles=100,
    x_entry=0,
    x_start=-50,
    max_steps=1000,
    temperature=293,
    One_Directional_Propagation=False,
    min_cos_theta=0.5,
    Fast_Ballistic_Propagation=True
):
   
    results = []
    for i in range(n_particles):
        
        p = Particle(
            name=particle.name,
            mass=particle.mass,
            charge=particle.charge,
            initial_velocity=list(particle.initial_velocity)
        )

        if simulation_mode == "thermal_only":
            
            res = propagate_particle(
                material, 
                p,
                x_entry=x_entry, 
                x_start=x_start,
                max_steps=max_steps,
                thermal_mode=True,
                One_Directional_Propagation=One_Directional_Propagation,
                min_cos_theta=min_cos_theta,
                Fast_Ballistic_Propagation=False
            )

        elif simulation_mode == "non_thermal_only":
            
            res = propagate_particle(
                material, 
                p,
                x_entry=x_entry, 
                x_start=x_start,
                max_steps=max_steps,
                thermal_mode=False,
                One_Directional_Propagation=One_Directional_Propagation,
                min_cos_theta=min_cos_theta,
                Fast_Ballistic_Propagation=True
            )

        elif simulation_mode == "full_moderation":
            res = propagate_moderated_particle(
                material, 
                p,
                x_entry=x_entry, 
                x_start=x_start,
                max_steps=max_steps,
                temperature=temperature,
                One_Directional_Propagation=One_Directional_Propagation,
                min_cos_theta=min_cos_theta,
                Fast_Ballistic_Propagation=Fast_Ballistic_Propagation
            )

        else:
            raise ValueError(f"Unknown simulation_mode: {simulation_mode!r}")

        results.append(res)

    return results




def get_material_for_two_slabs(x, 
                               mat1=Water, 
                               L1=Water.length_x, 
                               mat2=Graphite, 
                               L2=Graphite.length_x
                               ):



    if 0 <= x < L1:
        return mat1
    elif L1 <= x <= L1 + L2:
        return mat2
    else:
        return None
    
def propagate_two_slabs(
    mat1,
    mat2,
    L1,
    L2,
    particle,
    x_entry=0,
    x_start=-50,
    max_steps=1000,
    temperature=293,
    One_Directional_Propagation=False,
    min_cos_theta=0.5,
    Fast_Ballistic_Propagation=True
):
    
    pos = np.array([x_start, 0.0, 0.0])
    path = [pos.copy().tolist()]
    energies = [particle.energy_ev]
    thermalized_index = None
    energies_in_material = []
    end_position = None


    for i in range(1, 11):
        xi = x_start + (x_entry - x_start)*(i/10)
        pos = np.array([xi, 0.0, 0.0])
        path.append(pos.tolist())
        energies.append(particle.energy_ev)

    total_length = L1 + L2

    for _ in range(max_steps):
        x = pos[0]

        if 0 <= x < L1:
            material = mat1
        elif L1 <= x <= total_length:
            material = mat2
        else:
            print(f"Neutron exited slab at x = {pos[0]:.3f} m")
            end_position = x
            break


        if particle.is_thermal(temperature):
            if thermalized_index is None:
                thermalized_index = len(path) - 1
            direction = random_unit_vector(forward_only=One_Directional_Propagation)
            step = sample_step_length(material, energy_ev=None, fast=False)
        else:
            direction = forward_biased_unit_vector(min_cos_theta)
            step = sample_step_length(
                material,
                energy_ev=particle.energy_ev,
                fast=Fast_Ballistic_Propagation
            )


        pos = pos + step * direction
        path.append(pos.tolist())
        energies.append(particle.energy_ev)


        if 0 <= pos[0] <= total_length:
            energies_in_material.append(particle.energy_ev)


        p_abs = compute_absorption_probability(material, particle.energy_ev, step)
        if np.random.rand() < p_abs:
            print(f"Absorbed particle at x = {pos[0]:.3f} m, energy = {particle.energy_ev:.2f} eV")
            end_position = pos[0]
            break


        if not particle.is_thermal(temperature):
            particle.energy_ev *= material.energy_loss_factor
            particle.energy = particle.energy_ev * 1.602e-19

    path_with_energy = [p + [e] for p, e in zip(path, energies)]
    
    return path, path_with_energy, thermalized_index, energies, energies_in_material, end_position
