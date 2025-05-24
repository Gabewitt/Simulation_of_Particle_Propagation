from constants import Material, Particle



Water = Material(
    "Water",
    cross_section=169.0,          # m^-1
    density=1000.0,              # kg/m^3
    length_x=None,              # m
    absorption_probability=0.13,
    mass_number=18
)

Graphite = Material(
    "Graphite",
    cross_section=53,           # m^-1
    density=1700.0,              # kg/m^3
    length_x=None,              # m
    absorption_probability=0.0074,
    mass_number=12
)

Lead = Material(
    "Lead",
    cross_section=42.6,           # m^-1
    density=11340.0,             # kg/m^3
    length_x=None,              # m
    absorption_probability=0.13,
    mass_number=207
)

Air = Material(
    "Air",
    cross_section=0.03,          # m^-1
    density=1.2,                 # kg/m^3
    length_x=None,              # m
    absorption_probability=0.0001,
    mass_number=14.6
)

Borated_Polyethylene = Material(
    "Borated PE",
    cross_section=100,          # m^-1
    density=940.0,               # kg/m^3
    length_x=None,                # m
    absorption_probability=0.5,
    mass_number=10
)


# light water coolant for Nuclear reactors
LightWater = Water

# control rod (boron carbide) for Nuclear reactors
ControlRod = Material(
    "ControlRod",
    cross_section=384,    # m^-1
    density=2520.0,
    length_x=None,
    absorption_probability=0.99,
    mass_number=10
)


material_presets = {
    "Water": Water,
    "Graphite": Graphite,
    "Lead": Lead,
    "Air": Air,
    "Borated PE": Borated_Polyethylene,
    "LightWater": LightWater,
    "ControlRod": ControlRod
}




def get_predefined_material(name):
    return material_presets.get(name, None)



