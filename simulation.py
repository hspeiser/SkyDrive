import math
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------------
# This script simulates spinning a 2-meter long rod with the motor at its center.
# Then it calculates the total energy (and thus average power) consumed during spin-up
# by numerically integrating torque*omega over time.
# Finally, it simulates the rocket's vertical launch after release and plots the trajectory.
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
rod_diameter = 0.02    # assumed rod diameter in meters

Cd_rocket = 0.7
A_rocket = 0.01        # assumed rocket frontal area (m^2)

desired_tip_speed = 171.0  # m/s at radius 1m
spin_up_time = 300.0       # seconds (5 minutes)

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
# 3. Aerodynamic Drag Torque on the Rod
#    (We'll use a function so we can evaluate at any omega)
# ---------------------------------
N = 100
dr = radius / N
A_segment = rod_diameter * dr  # frontal area for each small segment along rod

def rod_drag_torque(omega):
    """
    Numerically integrate rod drag from r=0 to r=1, then double it.
    dτ = r * dF,  where dF = 0.5 * rho * Cd_rod * A_segment * (v^2)
    and v = omega*r for that segment.
    """
    torque_half = 0.0
    for i in range(N):
        r_mid = (i+0.5)*dr  # midpoint of the small segment
        v = omega * r_mid
        dF = 0.5 * rho * Cd_rod * A_segment * (v**2)
        dT = r_mid * dF
        torque_half += dT
    return 2.0 * torque_half  # Double for both sides of rod

rod_torque_drag_final = rod_drag_torque(omega_final)

print("--------------------------------------------------")
print("ROD AERODYNAMIC DRAG AT FINAL SPEED:")
print(f"Rod drag torque at final omega: {rod_torque_drag_final:.6f} N*m")

# ---------------------------------
# 4. Torque & Power at Final Speed (snapshot)
# ---------------------------------
alpha = omega_final / spin_up_time
T_inertia = I_total * alpha
T_final = T_inertia + rod_torque_drag_final  # if we were still accelerating at final speed
P_final = T_final * omega_final
P_maintain = rod_torque_drag_final * omega_final

print("--------------------------------------------------")
print("TORQUE AND POWER REQUIREMENTS (Snapshot at Final Speed):")
print(f"Spin-up time: {spin_up_time} s")
print(f"Angular acceleration (alpha): {alpha:.6f} rad/s^2")
print(f"Inertia torque (T_inertia) at final speed: {T_inertia:.6f} N*m")
print(f"Total torque at final speed (if still accelerating): {T_final:.6f} N*m")
print(f"Power if still accelerating at final speed: {P_final:.6f} W")
print(f"Maintenance torque at final speed (just drag): {rod_torque_drag_final:.6f} N*m")
print(f"Maintenance power at final speed (just drag): {P_maintain:.6f} W")

# ---------------------------------
# 5. *Total Energy* and *Average Power* Consumed During Spin-Up
#    (Numeric integration over time, assuming constant alpha)
# ---------------------------------

dt = 0.5  # time step for integration (can be smaller for more accuracy)
num_steps = int(spin_up_time / dt)
energy_consumed = 0.0

for i in range(num_steps+1):
    t = i * dt
    # Current angular velocity (assuming linear ramp-up from 0 to omega_final)
    omega_t = alpha * t
    # Current drag torque
    T_drag_t = rod_drag_torque(omega_t)
    # Inertial torque is constant if alpha is constant
    T_inertia_t = I_total * alpha
    # Total torque
    T_total_t = T_inertia_t + T_drag_t
    # Instantaneous power
    P_t = T_total_t * omega_t
    # Accumulate energy
    energy_consumed += P_t * dt

# Average power during spin-up:
average_power_spinup = energy_consumed / spin_up_time

# Compare with final rotational kinetic energy (just for reference):
rot_kinetic_energy_final = 0.5 * I_total * (omega_final**2)

energy_consumed_Wh = energy_consumed / 3600.0
average_power_spinup = energy_consumed / spin_up_time

# Compare with final rotational kinetic energy (just for reference):
rot_kinetic_energy_final = 0.5 * I_total * (omega_final**2)
rot_kinetic_energy_final_Wh = rot_kinetic_energy_final / 3600.0

print("--------------------------------------------------")
print("TOTAL ENERGY & AVERAGE POWER DURING SPIN-UP:")
print(f"Total spin-up energy consumed: {energy_consumed:,.2f} J "
      f"({energy_consumed_Wh:,.2f} Wh)")
print(f"Average power during spin-up: {average_power_spinup:,.2f} W")
print("")
print(f"Final rotational kinetic energy (for comparison): "
      f"{rot_kinetic_energy_final:,.2f} J "
      f"({rot_kinetic_energy_final_Wh:,.2f} Wh)")

# ---------------------------------
# 6. Simulate Rocket Ascent After Release
# ---------------------------------
time_step = 0.01
v = desired_tip_speed  # release speed (vertical)
h = 0.0
time_vals = []
height_vals = []
velocity_vals = []

t = 0.0
while v > 0:
    D = 0.5 * rho * Cd_rocket * A_rocket * v**2
    a = -g - (D / rocket_mass)
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
# 7. Plot the rocket trajectory
# ---------------------------------
plt.figure(figsize=(10,6))
plt.plot(time_vals, height_vals, label='Rocket Altitude')
plt.title("Rocket Altitude Over Time After Release")
plt.xlabel("Time (s)")
plt.ylabel("Altitude (m)")
plt.grid(True)
plt.legend()
plt.show()
