import math
import numpy as np
import matplotlib.pyplot as plt
import json

# Define material properties
materials = {
    'Aluminum': {'density': 2700, 'UTS': 300e6},
    'Steel': {'density': 7850, 'UTS': 400e6},
    'Titanium': {'density': 4500, 'UTS': 900e6},
    'Carbon Fiber': {'density': 1600, 'UTS': 600e6}
}

def read_input_from_json(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data

def main():
    # Read inputs from JSON file
    input_file = 'spinlaunch_input.json'  # Change this to your JSON file path
    inputs = read_input_from_json(input_file)

    # Extract inputs
    P_motor = inputs['motor_power']           # Motor power in watts
    GR = inputs['gear_ratio']                 # Gear ratio (dimensionless)
    P_vacuum = inputs['vacuum_strength']      # Vacuum strength in Pascals
    m_payload = inputs['payload_weight']      # Payload weight in kilograms
    Cd = inputs['coefficient_of_drag']        # Coefficient of drag of the payload
    A_p = inputs['payload_cross_section_area']# Payload cross-sectional area in square meters
    material_name = inputs['arm_material']    # Arm material name
    R = inputs['arm_length']                  # Arm length (radius) in meters
    r_arm = inputs['arm_cross_section_radius']# Arm cross-sectional radius in meters

    # Material properties
    if material_name not in materials:
        raise ValueError(f"Material '{material_name}' not recognized. Available materials: {list(materials.keys())}")
    material = materials[material_name]
    rho_material = material['density']
    UTS_material = material['UTS']

    # Compute arm mass
    A_arm = math.pi * r_arm**2  # Cross-sectional area of the arm
    L_arm = R                   # Length of the arm
    m_arm = rho_material * A_arm * L_arm

    # Compute moment of inertia
    I_arm = (1/3) * m_arm * R**2
    I_payload = m_payload * R**2
    I_total = I_arm + I_payload

    # Compute maximum angular velocity (ω_max) based on material strength
    omega_max = math.sqrt((UTS_material * A_arm) / (R * (m_payload + 0.5 * m_arm)))

    # Convert omega_max to RPM
    omega_max_rpm = omega_max * (60 / (2 * math.pi))

    # Compute time to reach ω_max based on motor power
    t_max = (I_total * omega_max**2) / (2 * P_motor)

    # Compute exit velocity
    v_exit = omega_max * R

    # Simulate projectile motion with air drag after release
    g = 9.81        # Acceleration due to gravity (m/s^2)
    rho_air = 1.225 # Air density at sea level (kg/m^3)

    # Initial conditions
    v = v_exit
    h = 0
    t = 0
    dt = 0.01  # Time step (seconds)

    # Lists to store trajectory data
    t_list = [t]
    h_list = [h]
    v_list = [v]

    print("\nCalculating trajectory...")

    # Simulate ascent until the payload starts descending
    while v > 0:
        # Air drag force
        F_drag = 0.5 * rho_air * v**2 * A_p * Cd

        # Net acceleration
        a = -g - (F_drag / m_payload)

        # Update velocity and position
        v += a * dt
        h += v * dt

        # Update time
        t += dt

        # Store data
        t_list.append(t)
        h_list.append(h)
        v_list.append(v)

    # Maximum height reached
    h_max = max(h_list)

    # Output results
    print("\n=== Simulation Results ===")
    print(f"Maximum speed inside (RPM): {omega_max_rpm:.2f} RPM")
    print(f"Exit velocity: {v_exit:.2f} m/s")
    print(f"Time to reach maximum speed: {t_max:.2f} seconds")
    print(f"Maximum height reached: {h_max:.2f} meters")

    # Plot trajectory
    plt.figure(figsize=(10, 6))
    plt.plot(t_list, h_list)
    plt.title('Payload Trajectory After Release')
    plt.xlabel('Time (s)')
    plt.ylabel('Height (m)')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
