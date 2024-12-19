import math
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------------
# This script simulates spinning a 2-meter long arm with the motor at its center, and a mass of 
# 0.5 kg (the rocket) on one end at radius = 1 m, and a 0.5 kg counterweight on the opposite end (also at 1 m).
#
# The goal:
# 1. Achieve a tip speed of Mach 0.5 (≈171 m/s) at a radius of 1 m.
# 2. Calculate the required RPM and angular velocity for that tip speed.
# 3. Estimate motor power requirements to spin up to that speed within 5 minutes (300 seconds), 
#    assuming direct drive and ignoring aerodynamic drag for simplicity.
# 4. Model the rocket's vertical ascent after release at 171 m/s (vertical), with a Cd = 0.7 
#    and plot how high it goes.
#
# Changes from previous scenario:
# - The arm is 2 meters total: 1 meter radius on each side of the motor.
# - Mass distribution:
#   * Rocket mass = 0.5 kg at 1 m radius on one side.
#   * Counterweight = 0.5 kg at 1 m radius on opposite side.
# - Thus total mass at radius = 1m is now 1.0 kg (0.5 kg rocket on one side, 0.5 kg counterweight on other).
#
# Assumptions and constants:
# - Air density: ρ = 1.225 kg/m^3
# - Gravity: g = 9.81 m/s²
# - Time to spin up: 300 s
# - Desired tip speed: 171 m/s
# - Radius: 1.0 m from center to mass
# - Cd of rocket = 0.7
# - Rocket reference area A = 0.01 m² (assumed)
#
# Steps:
# 1. Compute angular velocity (ω) and RPM needed for tip speed.
# 2. Compute moment of inertia with the masses on both ends and find required torque and power.
# 3. Simulate vertical ascent of rocket after release at 171 m/s, include drag and find max altitude.
# 4. Plot altitude over time.
#
# We'll print extensive logging and add comments throughout.
# ---------------------------------------------------------------------------------

# ---------------------------------
# Known constants and parameters
# ---------------------------------
g = 9.81               # gravitational acceleration (m/s^2)
rho = 1.225            # air density (kg/m^3)
rocket_mass = 0.5      # rocket mass (kg)
counter_mass = 0.5     # counterweight mass (kg)
Cd_rocket = 0.7        # Coefficient of drag for rocket
A_rocket = 0.01        # Assumed reference area (m^2)
desired_tip_speed = 400 # m/s
radius = 1.0            # meter (radius of the arm from center to mass)
spin_up_time = 300.0    # 5 minutes to reach target speed

# ---------------------------------
# 1. Compute required RPM for tip speed
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
# 2. Compute torque and power needed (ignoring drag)
# ---------------------------------
# Moment of inertia: For two equal masses, each at radius = 1 m:
# I = m*r² + m*r² = (0.5 kg)*(1 m²) + (0.5 kg)*(1 m²) = 1.0 kg*m² total
#
# Rod mass is negligible compared to these masses, so we ignore rod mass.
I = (rocket_mass * radius**2) + (counter_mass * radius**2)

print("--------------------------------------------------")
print("MOMENT OF INERTIA CALCULATION:")
print(f"Rocket mass: {rocket_mass:.3f} kg at {radius} m")
print(f"Counterweight mass: {counter_mass:.3f} kg at {radius} m")
print(f"Total I: {I:.6f} kg*m^2")

# Angular acceleration needed:
alpha = omega_final / spin_up_time  # rad/s²

# Torque needed (no drag):
tau = I * alpha

# Power at final speed (no drag):
power_final = tau * omega_final

print("--------------------------------------------------")
print("TORQUE AND POWER REQUIREMENTS (NO DRAG):")
print(f"Spin-up time: {spin_up_time} s")
print(f"Angular acceleration (alpha): {alpha:.6f} rad/s^2")
print(f"Torque (tau) required (no drag): {tau:.6f} N*m")
print(f"Power at final speed (no drag): {power_final:.6f} Watts")

# ---------------------------------
# 3. Simulate the rocket ascent after release:
# ---------------------------------
# The rocket is released vertically upwards at 171 m/s.
# We model vertical motion under gravity and drag:
#
# dv/dt = -g - (D/m), D = 0.5 * rho * Cd_rocket * A_rocket * v²
#
# We'll integrate until v <= 0 (rocket stops ascending).

time_step = 0.01  # small time step for simulation
v = desired_tip_speed
h = 0.0
time_vals = []
height_vals = []
velocity_vals = []

t = 0.0
while v > 0: 
    # Calculate drag (always opposes velocity)
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
print("Check the generated plot for altitude vs. time.")

# ---------------------------------
# 4. Plot the rocket trajectory
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
