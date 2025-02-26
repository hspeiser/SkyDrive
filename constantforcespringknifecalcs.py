import numpy as np
import matplotlib.pyplot as plt

def inch_to_meter(inch):
    """Convert inches to meters."""
    return inch * 0.0254

def rpm_to_rad_per_sec(rpm):
    """Convert revolutions per minute to radians per second."""
    return rpm * (2.0 * np.pi / 60.0)

def slender_rod_inertia(mass, length):
    """
    Moment of inertia of a slender rod of 'length' (pivot at one end)
    about that end. I = (1/3) * M * L^2.
    
    length here is the entire length of the rod in meters.
    If your pivot is not exactly at the end, you can split the rod into sections
    or use parallel-axis theorem. This is a simplified approach.
    """
    return (1.0/3.0) * mass * (length**2)

def required_torque(I, theta_radians, deployment_time):
    """
    Assuming constant angular acceleration from rest to the angle 'theta_radians' 
    in 'deployment_time', the angular acceleration alpha = theta / (t^2 / 2) * 2,
    or more simply alpha = 2 * theta / (deployment_time^2) if starting from 0 speed.
    
    Then the required *constant* torque Tau = I * alpha.
    """
    alpha = 2.0 * theta_radians / (deployment_time**2)
    Tau = I * alpha
    return Tau

def torque_to_force(torque, radius):
    """
    Given a torque and a radius (perpendicular distance from pivot 
    to the spring's line of action), compute the required force.
    
    Force = Torque / radius
    """
    # If radius is zero or extremely small, watch out for division by zero!
    if radius <= 0:
        return 0.0
    return torque / radius

def main():
    # ------------------------------
    # USER-DEFINED PARAMETERS
    # ------------------------------
    
    # Knife parameters
    knife_length_in = 7.5         # total length in inches
    pivot_from_one_end_in = 7.0   # pivot is 7" from that end
    mass_knife_kg = 0.2          # example mass of the knife (adjust to your actual mass)
    
    # Deployment parameters
    deployment_time_s = 0.05      # how quickly (in seconds) you want the knife to rotate 90°
    rotate_angle_deg = 90.0       # rotating from 0° to 90°
    
    # System spin info (if needed for context)
    # e.g. 1900 rpm about the central axis where the rope is 1 m away
    system_rpm = 1900
    system_omega_rad_s = rpm_to_rad_per_sec(system_rpm)
    
    # Convert units
    knife_length_m = inch_to_meter(knife_length_in)
    pivot_from_one_end_m = inch_to_meter(pivot_from_one_end_in)
    
    # For a slender rod pivot about one end:
    # If your pivot is truly at the end, you can use the standard I = (1/3) M L^2.
    # But here, the pivot is near one end, not exactly at the tip.
    # This code uses the simple formula for the entire length. 
    # Adapt the inertia if you want to be more precise.
    I_knife = slender_rod_inertia(mass_knife_kg, knife_length_m)
    
    # We want torque for 90° in 'deployment_time_s'.
    theta_radians = np.deg2rad(rotate_angle_deg)
    torque_req = required_torque(I_knife, theta_radians, deployment_time_s)
    
    # We'll sample distances in increments of 0.5 inch along the blade:
    # e.g., from 0.5" up to 7.5". (Adjust range as needed.)
    distances_in = np.arange(0.5, knife_length_in+0.1, 0.5)
    distances_m = inch_to_meter(distances_in)
    
    forces = []
    
    for r in distances_m:
        f = torque_to_force(torque_req, r)
        forces.append(f)
    
    # -------------
    # PRINT RESULTS
    # -------------
    print(f"System spin speed: {system_rpm} rpm = {system_omega_rad_s:.2f} rad/s (for reference)")
    print(f"Knife length = {knife_length_in} in = {knife_length_m:.4f} m")
    print(f"Pivot from end = {pivot_from_one_end_in} in = {pivot_from_one_end_m:.4f} m")
    print(f"Mass of knife = {mass_knife_kg} kg")
    print(f"Desired deployment: {rotate_angle_deg} deg in {deployment_time_s} s")
    print("---------------------------------------------------")
    print(f"Approx. moment of inertia (slender rod assumption): {I_knife:.6f} kg*m^2")
    print(f"Torque required for constant accel to 90 deg in {deployment_time_s}s: {torque_req:.6f} N*m")
    print("---------------------------------------------------")
    print("Distance from pivot (in) | Required spring force (N)")
    for (d_in, f_n) in zip(distances_in, forces):
        print(f"{d_in:>22.2f}          {f_n:>18.3f}")
    
    # -------------
    # PLOT RESULTS
    # -------------
    plt.figure(figsize=(8, 5))
    plt.plot(distances_in, forces, marker='o', linestyle='-')
    plt.title("Required Constant-Force Spring vs. Mount Radius")
    plt.xlabel("Mount distance from pivot (inches)")
    plt.ylabel("Required Force (N)")
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
