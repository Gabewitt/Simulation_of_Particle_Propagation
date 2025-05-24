import os
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go




def save_interactive_neutron_plot(
    path,
    material_length_x,
    filename="neutron_path.html",
    thermalized_index=None
):
    
    x, y, z = zip(*path)

    
    max_dev_y = max(abs(v) for v in y)
    max_dev_z = max(abs(v) for v in z)
    base_extent = max(max_dev_y, max_dev_z)
    wall_extent = base_extent * 1.1  
    if wall_extent < 1.0:
        wall_extent = 1.0

    # Create a meshgrid for the entry and exit walls
    Y = np.linspace(-wall_extent, wall_extent, 2)
    Z = np.linspace(-wall_extent, wall_extent, 2)
    Y_grid, Z_grid = np.meshgrid(Y, Z)
    X0 = np.zeros_like(Y_grid)
    X1 = np.full_like(Y_grid, material_length_x)

    fig = go.Figure()

    # Entry wall at x = 0
    fig.add_trace(go.Surface(
        x=X0, y=Y_grid, z=Z_grid,
        showscale=False,
        opacity=0.1,
        colorscale=[[0, 'red'], [1, 'red']],
        name="Entry Plane (x=0)"
    ))

    # Exit wall at x = material_length_x
    fig.add_trace(go.Surface(
        x=X1, y=Y_grid, z=Z_grid,
        showscale=False,
        opacity=0.1,
        colorscale=[[0, 'red'], [1, 'red']],
        name=f"Exit Plane (x={material_length_x})"
    ))

    # Neutron path
    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z,
        mode='lines+markers',
        marker=dict(size=2, color='cyan'),
        line=dict(color='blue', width=2),
        name="Neutron Path"
    ))

    # Start marker
    fig.add_trace(go.Scatter3d(
        x=[x[0]], y=[y[0]], z=[z[0]],
        mode='markers',
        marker=dict(size=6, color='green', symbol='circle'),
        name="Start"
    ))

    # End marker
    fig.add_trace(go.Scatter3d(
        x=[x[-1]], y=[y[-1]], z=[z[-1]],
        mode='markers',
        marker=dict(size=6, color='red', symbol='circle'),
        name="End"
    ))

    # Thermalization point, if any
    if thermalized_index is not None:
        fig.add_trace(go.Scatter3d(
            x=[x[thermalized_index]],
            y=[y[thermalized_index]],
            z=[z[thermalized_index]],
            mode='markers',
            marker=dict(size=6, color='blue', symbol='diamond'),
            name="Thermalization"
        ))

    # Center the scene around y=z=0 and give some margin on x
    fig.update_layout(
        scene=dict(
            xaxis=dict(title='x (m)', range=[-material_length_x*0.1, material_length_x*1.1]),
            yaxis=dict(title='y (m)', range=[-wall_extent, wall_extent]),
            zaxis=dict(title='z (m)', range=[-wall_extent, wall_extent]),
            camera=dict(
                center=dict(x=0.5, y=0, z=0),  # halfway between walls on x; y,z centered
                eye=dict(x=1.5, y=1.5, z=1.0)
            )
        ),
        title='Interactive Neutron Propagation',
        margin=dict(l=0, r=0, b=0, t=40),
        legend=dict(x=0.01, y=0.99)
    )

    fig.write_html(filename)
    
  

def save_interactive_neutron_plot_multiple(
    paths,                  
    material_length_x,
    filename="multiple_neutron_simulation.html",
    thermalized_indices=None
):
   
    all_y = np.concatenate([np.array([p[1] for p in path]) for path in paths])
    all_z = np.concatenate([np.array([p[2] for p in path]) for path in paths])
    base_extent = max(all_y.max(), -all_y.min(), all_z.max(), -all_z.min())
    wall_extent = max(base_extent * 1.1, 1.0)

    
    Y = np.linspace(-wall_extent, wall_extent, 2)
    Z = np.linspace(-wall_extent, wall_extent, 2)
    Yg, Zg = np.meshgrid(Y, Z)
    X0 = np.zeros_like(Yg)
    X1 = np.full_like(Yg, material_length_x)

    fig = go.Figure()

    
    fig.add_trace(go.Surface(x=X0, y=Yg, z=Zg,
                             showscale=False, opacity=0.1,
                             colorscale=[[0,'red'],[1,'red']],
                             name="Entry (x=0)"))
    fig.add_trace(go.Surface(x=X1, y=Yg, z=Zg,
                             showscale=False, opacity=0.1,
                             colorscale=[[0,'red'],[1,'red']],
                             name=f"Exit (x={material_length_x})"))

    
    for i, path in enumerate(paths):
        x, y, z = zip(*path)
        fig.add_trace(go.Scatter3d(
            x=x, y=y, z=z,
            mode='lines',
            line=dict(width=1),
            name=f"Neutron {i+1}"
        ))
        
        if thermalized_indices and thermalized_indices[i] is not None:
            idx = thermalized_indices[i]
            fig.add_trace(go.Scatter3d(
                x=[x[idx]], y=[y[idx]], z=[z[idx]],
                mode='markers',
                marker=dict(size=4, symbol='diamond', color='blue'),
                name=f"Therm {i+1}"
            ))

    
    fig.update_layout(
        scene=dict(
            xaxis=dict(title='x (m)', range=[-material_length_x*0.1, material_length_x*1.1]),
            yaxis=dict(title='y (m)', range=[-wall_extent, wall_extent]),
            zaxis=dict(title='z (m)', range=[-wall_extent, wall_extent]),
            camera=dict(center=dict(x=0.5, y=0, z=0),
                        eye=dict(x=1.5, y=1.5, z=1.0))
        ),
        title='Multiple Neutron Propagation',
        margin=dict(l=0, r=0, b=0, t=40)
    )

    fig.write_html(filename)


def save_abs_trans_fraction_plot(results, output_folder, base_filename):
    
    n = len(results)
    
    lengths = [len(r[4]) for r in results]
    max_len = max(lengths)

    steps = list(range(1, max_len+1))
    frac_transmitted = []
    for i in range(max_len):
        
        count_trans = sum(1 for ln in lengths if ln > i)
        frac_transmitted.append(count_trans / n)

    frac_absorbed = [1 - ft for ft in frac_transmitted]

    
    plt.figure(figsize=(8,5), dpi=500)
    plt.plot(steps, frac_absorbed, label='Fraction Absorbed')
    plt.plot(steps, frac_transmitted, label='Fraction Transmitted')
    plt.xlabel("Step Index")
    plt.ylabel("Fraction of Neutrons")
    plt.title("Absorption vs Transmission over Steps")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    fname = os.path.join(output_folder, f"{base_filename}_fraction.png")
    plt.savefig(fname)
    plt.close()
    return fname





def save_interactive_two_slab_plot(
    path,
    L1,
    L2,
    filename="neutron_two_slab.html",
    thermalized_index=None
):
    
    x, y, z = zip(*path)

    
    max_dev_y = max(abs(v) for v in y)
    max_dev_z = max(abs(v) for v in z)
    wall_extent = max(max_dev_y, max_dev_z) * 1.1
    wall_extent = max(wall_extent, 1.0)

    
    Y = np.linspace(-wall_extent, wall_extent, 2)
    Z = np.linspace(-wall_extent, wall_extent, 2)
    Yg, Zg = np.meshgrid(Y, Z)

    
    planes = [
        (0.0, "Entry (x=0)", "red"),
        (L1,   f"Interface (x={L1})", "orange"),
        (L1+L2, f"Exit (x={L1+L2})", "red")
    ]

    fig = go.Figure()

    
    for xpos, name, color in planes:
        Xg = np.full_like(Yg, xpos)
        fig.add_trace(go.Surface(
            x=Xg, y=Yg, z=Zg,
            showscale=False,
            opacity=0.2,
            colorscale=[[0, color],[1, color]],
            name=name,
            showlegend=True
        ))

    
    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z,
        mode='lines+markers',
        marker=dict(size=2, color='cyan'),
        line=dict(width=2, color='blue'),
        name="Neutron Path"
    ))

    
    fig.add_trace(go.Scatter3d(
        x=[x[0]], y=[y[0]], z=[z[0]],
        mode='markers', marker=dict(size=6, color='green'),
        name="Start"
    ))
    fig.add_trace(go.Scatter3d(
        x=[x[-1]], y=[y[-1]], z=[z[-1]],
        mode='markers', marker=dict(size=6, color='red'),
        name="End"
    ))

    
    if thermalized_index is not None:
        fig.add_trace(go.Scatter3d(
            x=[x[thermalized_index]],
            y=[y[thermalized_index]],
            z=[z[thermalized_index]],
            mode='markers',
            marker=dict(size=6, color='blue', symbol='diamond'),
            name="Thermalized"
        ))

   
    total_length = L1 + L2
    fig.update_layout(
        scene=dict(
            xaxis=dict(title='x', range=[-0.1*total_length, 1.1*total_length]),
            yaxis=dict(title='y', range=[-wall_extent, wall_extent]),
            zaxis=dict(title='z', range=[-wall_extent, wall_extent]),
            camera=dict(eye=dict(x=1.5, y=1.2, z=1.0))
        ),
        title="Neutron Through Two Slabs",
        margin=dict(l=0, r=0, b=0, t=30)
    )

    fig.write_html(filename)
    return filename