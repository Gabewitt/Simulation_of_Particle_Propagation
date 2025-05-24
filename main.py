import numpy as np
import os
from datetime import datetime

from constants import Material, Particle


from properties import (
Water, Graphite, Lead, Air,
Borated_Polyethylene, LightWater, ControlRod,
get_predefined_material
)



from outputs import (
load_config_yaml,
save_path_to_csv,
save_multiple_paths_to_csv,
write_summary_file
)


from systems import (
propagate_particle,
propagate_moderated_particle,
simulate_multiple_neutrons,
propagate_two_slabs
)


from visualisations import (
save_interactive_neutron_plot,
save_interactive_neutron_plot_multiple, save_abs_trans_fraction_plot,
save_interactive_two_slab_plot
)


def main():
    
    # Load configuration from YAML file
    
    cfg = load_config_yaml("simulation_config.yaml")
    
    if cfg["material_choice"].lower() == "custom":
        m = cfg["custom_material"]
        material = Material(
            name=m["name"],
            cross_section=m["cross_section"],
            density=m["density"],
            length_x=m["length_x"],
            absorption_probability=m["absorption_probability"],
            mass_number=m["mass_number"]
        )
    else:
        material = get_predefined_material(cfg["material_choice"])
        if material is None:
            raise ValueError(f"Material '{cfg['material_choice']}' not found in preperties.py.")
        material.length_x = cfg["material_length"]
    
    if cfg["simulation_two_material"]:
        material_1 = get_predefined_material(cfg["material_choice_1"])
        material_2 = get_predefined_material(cfg["material_choice_2"])
        
        if material_1 is None or material_2 is None:
            raise ValueError(f"Material '{cfg['material_choice_1']}' or '{cfg['material_choice_2']}' not found in preperties.py.")
        material_1.length_x = cfg["material_length_1"]
        material_2.length_x = cfg["material_length_2"]
        
    
    
    p = cfg["particle"]
    initial_velocity = [float(v) for v in p["initial_velocity"]]
    particle = Particle(
        name=p["name"],
        mass=p["mass"],
        charge=p["charge"],
        initial_velocity=initial_velocity
    )

    now = datetime.now().strftime("%Y-%m-%d_%H-%M")
    output_root = cfg["output"].get("output_dir", "Propagation_Data")
    os.makedirs(output_root, exist_ok=True)
    base_filename = f"{material.name}_{particle.name}_{cfg['simulation_mode']}_{now}"
    output_folder = os.path.join(output_root, base_filename)
    os.makedirs(output_folder, exist_ok=True)
    
    sim = cfg["simulation"]
    sim["x_start"] = float(sim["x_start"])
    sim["x_entry"] = float(sim["x_entry"])
    flags = cfg["flags"]
    mode = cfg["simulation_mode"]
    number_of_simulations = sim["number_of_simulations"]
    
    
    
    if number_of_simulations == 1: # Single particle simulation
        
    
        if mode == "thermal_only":
            
            particle.set_thermal_velocity(sim["temperature"])
            
            path, path_with_energy, thermalized_index, _, _, _ = propagate_particle(
                material,
                particle,
                x_entry=sim["x_entry"],
                x_start=sim["x_start"],
                max_steps=sim["max_steps"],
                thermal_mode=True,
                One_Directional_Propagation=flags["One_Directional_Propagation"],
                min_cos_theta=flags["min_cos_theta"],
                Fast_Ballistic_Propagation=False
            )
            
            

        
        elif mode == "non_thermal_only":
            path, path_with_energy, thermalized_index, _, _, _ = propagate_particle(
                material,
                particle,
                x_entry=sim["x_entry"],
                x_start=sim["x_start"],
                max_steps=sim["max_steps"],
                thermal_mode=False,
                One_Directional_Propagation=flags["One_Directional_Propagation"],
                min_cos_theta=flags["min_cos_theta"],
                Fast_Ballistic_Propagation=True
            )
            
        
        elif mode == "full_moderation":
            
            if cfg["simulation_two_material"]:
    
                path, path_with_energy, thermal_idx, _, _, _ = propagate_two_slabs(
                    mat1=material_1, 
                    mat2=material_2, 
                    L1=material_1.length_x, 
                    L2=material_2.length_x,
                    particle=particle, 
                    x_entry=sim["x_entry"], 
                    x_start=sim["x_start"],
                    max_steps=sim["max_steps"], 
                    temperature=sim["temperature"],
                    One_Directional_Propagation=flags["One_Directional_Propagation"],
                    min_cos_theta=flags["min_cos_theta"],
                    Fast_Ballistic_Propagation=flags["Fast_Ballistic_Propagation"]
                )
            else:
                
                path, path_with_energy, thermalized_index, _, _, _ =propagate_moderated_particle(
                    material,
                    particle,
                    x_entry=sim["x_entry"],
                    x_start=sim["x_start"],
                    max_steps=sim["max_steps"],
                    temperature=sim["temperature"],
                    One_Directional_Propagation=flags["One_Directional_Propagation"],
                    min_cos_theta=flags["min_cos_theta"],
                    Fast_Ballistic_Propagation=flags["Fast_Ballistic_Propagation"]
                )
            
        else:
            raise ValueError(f"Unknown simulation mode: {mode}")


        # Save results
        
        if cfg["simulation_two_material"]:
            csv_path     = os.path.join(output_folder, f"{base_filename}_2slabs.csv")
            html_path    = os.path.join(output_folder, f"{base_filename}_2slabs.html")
            summary_path = os.path.join(output_folder, f"{base_filename}_2slabs.txt")

            if cfg["output"]["save_csv"]:
                save_path_to_csv(path_with_energy, csv_path)
            if cfg["output"]["save_plot"]:
                save_interactive_two_slab_plot(
                    path, L1=material_1.length_x, L2=material_2.length_x,
                    filename=html_path,
                    thermalized_index=thermal_idx
                )

            write_summary_file(
                summary_path,
                [material_1, material_2],
                particle,
                cfg,
                html_path=html_path,
                csv_path=csv_path,
                date_time=now,
                path=path,
                thermalized_index=thermal_idx
            )

            print(f"Two-slab run complete.  CSV to {csv_path}, HTML to {html_path}, TXT to {summary_path}")
        else:
            csv_path = os.path.join(output_folder, f"{base_filename}.csv")
            html_path = os.path.join(output_folder, f"{base_filename}.html")
            summary_path = os.path.join(output_folder, f"{base_filename}.txt")

            
            
            if cfg["output"]["save_csv"]:
                save_path_to_csv(path_with_energy, csv_path)
            
            if cfg["output"]["save_plot"]:
                save_interactive_neutron_plot(path, material_length_x= material.length_x, filename=html_path, thermalized_index=thermalized_index)
                
            
            write_summary_file(summary_path, material, particle, cfg, html_path=html_path, csv_path= csv_path , date_time=now, path=path, thermalized_index=thermalized_index)
                

            
            
            print(f"Simulation '{mode}' with {material.name}  and with {particle.name} complete.")
    
    
    # Multiple particle simulation
    
    elif number_of_simulations > 1:
        results = simulate_multiple_neutrons(
            material,
            particle=particle,
            simulation_mode=mode,
            n_particles=number_of_simulations,
            x_entry=sim["x_entry"],
            x_start=sim["x_start"],
            max_steps=sim["max_steps"],
            temperature=sim["temperature"],
            One_Directional_Propagation=flags["One_Directional_Propagation"],
            min_cos_theta=flags["min_cos_theta"],
            Fast_Ballistic_Propagation=flags["Fast_Ballistic_Propagation"]
        )
     
        # Save results for multiple particles
    
        paths               = [res[0] for res in results]
        paths_with_energy   = [res[1] for res in results]  
        thermalized_indices = [res[2] for res in results]  
        csv_multiple     = os.path.join(output_folder, f"{base_filename}_multiple.csv")
        html_multiple    = os.path.join(output_folder, f"{base_filename}_multiple.html")
        summary_multiple = os.path.join(output_folder, f"{base_filename}_multiple.txt")
        
        if cfg["output"]["save_csv"]:
            save_multiple_paths_to_csv(paths_with_energy, csv_multiple)
        if cfg["output"]["save_plot"]:
            save_interactive_neutron_plot_multiple(
                paths, 
                material_length_x= material.length_x, filename=html_multiple, thermalized_indices=thermalized_indices)
        
        write_summary_file(
            summary_multiple,
            material,
            particle,
            cfg,
            html_path=html_multiple,
            csv_path=csv_multiple,
            date_time=now,
            path=None,               
            thermalized_index=None    
        )
        if cfg["output"]["save_absorption_transmission"]:
            fraction_png = save_abs_trans_fraction_plot(
                results,
                output_folder,
                base_filename
            )
            print(f"Saved absorption/transmission plot to {fraction_png}")
    
    else:
        raise ValueError(f"number_of_simulations must be a positive integer number.")
        
        
   

if __name__ == "__main__":
    main()
    
    
    
    
    
    
    
    