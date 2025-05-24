
import os
import csv
import yaml
import datetime
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import numpy as np


def load_config_yaml(filename):
    with open(filename, 'r') as f:
        config = yaml.safe_load(f)
    return config



def save_path_to_csv(path, filename):
    
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["x", "y", "z", "energy_eV"])
        writer.writerows(path)
    return filename


def save_multiple_paths_to_csv(paths, filename):
    
    n = len(paths)
    lengths = [len(p) for p in paths]
    max_len = max(lengths)

    headers = []
    for i in range(n):
        headers += [f"x_{i+1}", f"y_{i+1}", f"z_{i+1}", f"E_{i+1}_eV"]

    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)

        for step in range(max_len):
            row = []
            for p in paths:
                if step < len(p):
                    row += p[step]
                else:
                    row += ["", "", "", ""]
            writer.writerow(row)
            
            
def write_summary_file(
    filename,
    material,
    particle,
    cfg,
    html_path,
    csv_path,
    date_time=None,
    path=None,
    thermalized_index=None
):
    with open(filename, 'w') as f:
        f.write(f"Simulation Summary - {date_time}\n\n")

        f.write("Materials:\n")
        if isinstance(material, list):
            
            for m in material:
                f.write(f"  Name: {m.name}\n")
                f.write(f"    Cross Section: {m.cross_section} m^-1\n")
                f.write(f"    Density: {m.density} kg/m^3\n")
                f.write(f"    Absorption Probability: {m.absorption_probability}\n")
                f.write(f"    Mass Number: {m.mass_number}\n\n")
        else:
            f.write(f"  Name: {material.name}\n")
            f.write(f"  Cross Section: {material.cross_section} m^-1\n")
            f.write(f"  Density: {material.density} kg/m^3\n")
            f.write(f"  Absorption Probability: {material.absorption_probability}\n")
            f.write(f"  Mass Number: {material.mass_number}\n\n")

        f.write("Particle:\n")
        f.write(f"  Mass: {particle.mass} kg\n")
        f.write(f"  Charge: {particle.charge}\n")
        f.write(f"  Initial Velocity: {particle.initial_velocity} m/s\n")
        f.write(f"  Initial Energy: {particle.energy_ev:.4f} eV\n\n")

        sim = cfg["simulation"]
        f.write("Simulation Setup:\n")
        f.write(f"  Mode: {cfg['simulation_mode']}\n")
        f.write(f"  Temperature: {sim['temperature']} K\n")
        f.write(f"  Max Steps: {sim['max_steps']}\n")
        f.write(f"  x_start: {sim['x_start']} m\n")
        f.write(f"  x_entry: {sim['x_entry']} m\n")
        f.write(f"Number of Simulations: {sim['number_of_simulations']}\n\n")

        flags = cfg["flags"]
        f.write("Flags:\n")
        for key, val in flags.items():
            f.write(f"  {key}: {val}\n")

        f.write("\nOutputs:\n")
        if cfg["output"]["save_csv"]:
            f.write(f"  CSV File: {csv_path}\n")
        if cfg["output"]["save_plot"]:
            f.write(f"  Interactive HTML: {html_path}\n")

        
        if path is not None and thermalized_index is not None:
            x_th = path[thermalized_index][0]
            f.write(f"\nThermalization: neutron reached thermal energy at x = {x_th:.3f} m.\n")
        else:
            f.write("\nThermalization: neutron did NOT thermalize within the simulation domain.\n")

        if path is not None:
            x_final = path[-1][0]

            
            if isinstance(material, list):
                total_length = sum(m.length_x for m in material)
            else:
                total_length = material.length_x

            if 0 <= x_final <= total_length:
                f.write(f"Absorption: neutron was absorbed at x = {x_final:.3f} m.\n")
            else:
                f.write("Absorption: neutron exited the material without being absorbed.\n")