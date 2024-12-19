import math
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------------
# This script simulates spinning a 2-meter long rod with the motor at its center.
# The setup:
# - A 2 meter rod (1 meter radius on each side).
# - At one end (radius = 1 m), there is a 0.5 kg rocket.
# - At the other end (radius = 1 m), there is a 0.5 kg counterweight.
# - The rod itself weighs 0.944 kg (assume uniform).
#
# Goal:
# 1. Achieve a tip speed of Mach 0.5 (~171 m/s) at a radius of 1 m.
# 2. Consider the aerodynamic drag of the rod itself (Cd = 0.1), 
#    and include the torque required to overcome this drag.
#    We'll assume a certain rod diameter to calculate frontal area. Since not specified,
#    we must assume something reasonable. Let's assume a cylindrical rod of diameter 0.02 m.
#
# 3. Estimate the total torque and power needed by the motor to both accelerate the system to 
#    the final speed in under 5 minutes (300 s) AND overcome rod aerodynamic drag at the final speed.
#
# 4. Once released at 171 m/s vertically, simulate the rocket ascent under gravity and drag (Cd=0.7)
#    and plot how high it goes.
#
# Assumptions and approach:
# - Air density: ρ = 1.225 kg/m^3
# - Gravity: g = 9.81 m/s²
# - Spin-up time: 300 s
# - Desired tip speed = 171 m/s
# - Radius = 1 m
#
# Moment of Inertia:
# - Masses at the ends: rocket (0.5 kg) and counterweight (0.5 kg) both at 1 m radius.
#   I_masses = (0.5 * 1²) + (0.5 * 1²) = 1.0 kg*m²
#
# - Rod mass distribution:
#   The rod is 2 m long, mass = 0.944 kg, uniform. Moment of inertia about center:
#   I_rod = (1/12)*M*L² = (1/12)*0.944*(2²) = (1/12)*0.944*4 = (0.944*4)/12 = 3.776/12 = 0.314666... kg*m²
#
# Total I = I_masses + I_rod
#
# Aerodynamic drag on rod:
# - Assume rod is perpendicular to the direction of motion at all points (since it's spinning).
# - We'll integrate drag along the rod length. Each small segment at radius r sees 
#   v(r) = ω * r, and drag force dF = 0.5 * ρ * Cd_rod * A_segment * v(r)².
# - A_segment = rod_diameter * segment_length
# - Then torque contribution from that segment dτ = r * dF.
#
# We do a numeric integration along half the rod (0 to 1 m) and double it (since symmetrical).
#
# Finally, we find the torque required to maintain final angular velocity against drag (at final speed).
# The motor must provide at least this torque continuously at final speed, plus the torque to accelerate 
# up to that speed within 5 minutes.
#
# We'll print out all details and plot rocket trajectory.
#
# ---------------------------------------------------------------------------------

# ---------------------------------
# Known constants and parameters
# ---------------------------------
g = 9.81
rho = 1.225

rocket_mass = 0.5      # kg
counter_mass = 0.5     # kg
rod_mass = 0.944       # kg
Cd_rod = 0.1           # rod drag coefficient
rod_length = 2.0       # meters total
radius = 1.0           # meters on each side from center
rod_diameter = 0.02    # assumed rod diameter in meters for frontal area calc

Cd_rocket = 0.7
A_rocket = 0.01         # assumed rocket frontal area (m²)
desired_tip_speed = 200 # m/s at radius 1m
spin_up_time = 300.0     # s (5 minutes)

# ---------------------------------
# 1. Compute required angular velocity and RPM
# ---------------------------------
omega_final = desired_tip_speed / radius   # rad/s
rpm_final = (omega_final / (2*math.pi)) * 60

print("--------------------------------------------------")
print("COMPUTATION OF REQUIRED RPM:")
print(f"Desired tip speed: {desired_tip_speed} m/s")
print(f"Radius: {radius} m")
print(f"Angular velocity (omega_final): {omega_final:.3f} rad/s")
print(f"Required RPM: {rpm_final:.3f} RPM")

# ---------------------------------
# 2. Moment of Inertia
# ---------------------------------
# Masses at ends:
I_masses = (rocket_mass * radius**2) + (counter_mass * radius**2)
# Rod:
I_rod = (1.0/12.0)*rod_mass*(rod_length**2)
I_total = I_masses + I_rod

print("--------------------------------------------------")
print("MOMENT OF INERTIA:")
print(f"Mass at ends: rocket = {rocket_mass} kg, counter = {counter_mass} kg at 1 m")
print(f"I_masses = {I_masses:.6f} kg*m^2")
print(f"Rod mass = {rod_mass} kg, length = {rod_length} m")
print(f"I_rod = {I_rod:.6f} kg*m^2")
print(f"Total I = {I_total:.6f} kg*m^2")

# ---------------------------------
# 3. Rod Aerodynamic Drag Torque at final speed
# ---------------------------------
# We integrate along half the rod from 0 to 1m, then double.
# dτ = r * dF, dF = 0.5 * rho * Cd_rod * A_segment * (omega_final*r)^2
# A_segment = rod_diameter * segment_length
# We'll use numerical integration with N segments:
N = 100
dr = radius / N
A_segment = rod_diameter * dr

def rod_drag_torque(omega):
    torque_half = 0.0
    for i in range(N):
        r_mid = (i+0.5)*dr  # midpoint of segment
        v = omega * r_mid
        dF = 0.5 * rho * Cd_rod * A_segment * (v**2)
        dT = r_mid * dF
        torque_half += dT
    return 2.0 * torque_half  # double for both sides of rod

rod_torque_drag_final = rod_drag_torque(omega_final)

print("--------------------------------------------------")
print("ROD AERODYNAMIC DRAG AT FINAL SPEED:")
print(f"Rod drag torque at final omega: {rod_torque_drag_final:.6f} N*m")

# ---------------------------------
# 4. Torque and Power Requirements
# ---------------------------------
# Angular acceleration needed:
alpha = omega_final / spin_up_time

# In reality, torque must overcome inertia and drag that grows with speed.
# The final torque needed at steady speed (no acceleration) is just drag torque.
# But to accelerate, initially drag is negligible at low speed, grows with v².
# To find a rough estimate of required power:
#
# At final speed, to just maintain speed (no acceleration): 
#   T_maintain = rod_torque_drag_final  (since rocket and counterweight are points, negligible frontal area?)
#   P_maintain = T_maintain * omega_final
#
# To accelerate to final speed in spin_up_time:
#   T_inertia = I_total * alpha
#
# The total torque at final speed if we were still accelerating would be:
#   T_final_accel = T_inertia + rod_torque_drag_final
#
# However, after reaching final speed, acceleration stops, so the motor only needs drag torque to maintain.
# The maximum power output during acceleration will be somewhere before final speed, but let's 
# at least print the final values.

T_inertia = I_total * alpha
T_final = T_inertia + rod_torque_drag_final
P_final = T_final * omega_final  # if we were still accelerating at final speed

# Also compute just maintenance power at final speed (no acceleration):
P_maintain = rod_torque_drag_final * omega_final

print("--------------------------------------------------")
print("TORQUE AND POWER REQUIREMENTS:")
print(f"Spin-up time: {spin_up_time} s")
print(f"Angular acceleration (alpha): {alpha:.6f} rad/s^2")
print(f"Inertia torque (T_inertia) to accelerate: {T_inertia:.6f} N*m")
print(f"Total torque at final speed if still accelerating: {T_final:.6f} N*m")
print(f"Power if still accelerating at final speed: {P_final:.6f} W")
print(f"Maintenance torque at final speed (just drag): {rod_torque_drag_final:.6f} N*m")
print(f"Maintenance power at final speed (just drag): {P_maintain:.6f} W")

# Note: Real scenario requires integrating over speed to find exact energy and average power. 
# Here we provide a rough snapshot.

# ---------------------------------
# 5. Simulate Rocket Ascent After Release
# ---------------------------------
# Released at 171 m/s vertically, with Cd_rocket=0.7 and mass=0.5 kg:
#
# dv/dt = -g - (D/m), D = 0.5 * rho * Cd_rocket * A_rocket * v²
# Integrate until v ≤ 0.
time_step = 0.01
v = desired_tip_speed
h = 0.0
time_vals = []
height_vals = []
velocity_vals = []

t = 0.0
while v > 0: 
    D = 0.5 * rho * Cd_rocket * A_rocket * v**2
    a = -g - (D/rocket_mass)
    v = v + a*time_step
    h = h + v*time_step
    time_vals.append(t)
    height_vals.append(h)
    velocity_vals.append(v)
    t += time_step

max_height = max(height_vals)

print("--------------------------------------------------")
print("ROCKET TRAJECTORY AFTER RELEASE:")
print(f"Max height reached by rocket: {max_height:.3f} m")

# ---------------------------------
# 6. Plot the rocket trajectory
# ---------------------------------
plt.figure(figsize=(10,6))
plt.plot(time_vals, height_vals, label='Rocket Altitude')
plt.title("Rocket Altitude Over Time After Release")
plt.xlabel("Time (s)")
plt.ylabel("Altitude (m)")
plt.grid(True)
plt.legend()
plt.show()

# Done.
# ---------------------------------------------------------------------------------
