import math
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------------
# SPINNING "AIRFOIL" ROD SCRIPT
#   1) Prints out system details (mass, inertia, motor specs).
#   2) Uses a simpler airfoil-like drag model with a chord-based integration.
#   3) For each gear ratio, finds equilibrium angular velocity from:
#        T_motor_out = T_drag(omega)
#      Then checks feasibility vs. motor's no-load speed.
#   4) Prints everything in a table, and plots (a) Final RPM vs. Ratio,
#      and (b) Drag torque vs. omega for reference.
# ---------------------------------------------------------------------------------

# ------------------------
# 1) Setup / Constants
# ------------------------
rho = 1.225  # air density, kg/m^3
g = 9.81     # gravity, m/s^2

# System geometry
rod_length_total = 2.0   # meters, total (1 m each side)
radius = 1.0             # m from center to tip
chord = 0.02             # "airfoil" chord (m)
Cd_rod = 0.1             # drag coefficient for the rod as an airfoil

# Masses
rocket_mass = 0.5   # kg on one end
counter_mass = 0.5  # kg on the other end
rod_mass = 0.167    # kg (uniform rod)

# Motor specs
I_max = 200.0      # A
kv = 560.0         # RPM/V
Kt = 9.5493 / kv   # Nm/A
T_motor_max = I_max * Kt

# Battery/voltage assumption for motor no-load speed
#   For example:  12S @ 3.7 nominal => ~44.4 V
#   But let's just pick some voltage. 
#   We'll assume 44.4 V here for demonstration:
V_batt = 44.4
motor_no_load_rpm = kv * V_batt       # in RPM
motor_no_load_omega = (motor_no_load_rpm * 2.0 * math.pi) / 60.0  # rad/s

# ------------------------
# 2) Moment of Inertia
# ------------------------
I_masses = (rocket_mass * radius**2) + (counter_mass * radius**2)
I_rod = (1.0/12.0)*rod_mass*(rod_length_total**2)
I_total = I_masses + I_rod

# ------------------------
# 3) Print system details
# ------------------------
print("--------------------------------------------------")
print("SYSTEM SETUP & PARAMETERS (AIRFOIL DRAG MODEL):")
print(f"Rocket Mass:          {rocket_mass} kg")
print(f"Counterweight Mass:   {counter_mass} kg")
print(f"Rod Mass:             {rod_mass} kg")
print(f"Rod total length:     {rod_length_total} m (1 m each side)")
print(f"Airfoil chord:        {chord:.3f} m")
print(f"Airfoil Cd:           {Cd_rod}")
print("--------------------------------------------------")

print("MOMENT OF INERTIA:")
print(f"I_masses (ends):      {I_masses:.6f} kg*m^2")
print(f"I_rod (uniform):      {I_rod:.6f} kg*m^2")
print(f"I_total:              {I_total:.6f} kg*m^2")
print("--------------------------------------------------")

print("MOTOR SPECS:")
print(f"Motor kv:             {kv:.1f} RPM/V")
print(f"Battery Voltage:      {V_batt:.1f} V")
print(f"No-load motor RPM:    {motor_no_load_rpm:.1f} RPM")
print(f"No-load motor omega:  {motor_no_load_omega:.1f} rad/s")
print(f"I_max:                {I_max:.1f} A")
print(f"Kt (Nm/A):            {Kt:.6f} ")
print(f"Max motor torque:     {T_motor_max:.6f} Nm")
print("--------------------------------------------------")

# ------------------------
# 4) Airfoil-like Drag Integration
# ------------------------
# We'll integrate from r=0 to r=1 for half the rod, then double for both sides.
#   dA = chord * dr
#   dF = 0.5 * rho * Cd_rod * (v^2) * dA
#   v = omega * r
#   dT = r * dF
# We'll do numeric integration.

def rod_drag_torque_airfoil(omega):
    N = 200
    dr = radius / N
    torque_half = 0.0
    for i in range(N):
        r_mid = (i + 0.5)*dr
        v = omega * r_mid
        # local area
        dA = chord * dr
        dF = 0.5 * rho * Cd_rod * (v**2) * dA
        dT = r_mid * dF
        torque_half += dT
    return 2.0 * torque_half  # double for the other half

# Quick check of drag torque at some speeds
omegas_check = [0, 10, 50, 100, 200]
print("DRAG TORQUE CHECK (AIRFOIL MODEL):")
for omg in omegas_check:
    T_drag = rod_drag_torque_airfoil(omg)
    print(f"omega = {omg:3d} rad/s => Drag torque ~ {T_drag:.6f} N*m")
print("--------------------------------------------------")


# ------------------------
# 5) Solve T_motor_out = T_drag(omega) for each gear ratio
#    then check feasibility vs. no-load speed
# ------------------------
def find_equilibrium_omega_airfoil(T_req, omega_max=2000.0, tol=0.01):
    """
    Numerically solves rod_drag_torque_airfoil(omega) = T_req
    via a simple binary search between [0, omega_max].
    Returns the found 'omega' in rad/s, or clamp at omega_max if needed.
    """
    lower = 0.0
    upper = omega_max
    
    # Check boundary:
    if rod_drag_torque_airfoil(upper) < T_req:
        # Even at upper, we can't match T_req => can spin "faster" than we tested
        return upper
    
    while (upper - lower) > tol:
        mid = 0.5*(upper + lower)
        T_drag_mid = rod_drag_torque_airfoil(mid)
        if T_drag_mid < T_req:
            lower = mid
        else:
            upper = mid
    return 0.5*(upper + lower)

gear_ratio_list = [1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20]  # some plausible & not-so-plausible

equilibrium_rpm_list = []
feasible_list = []

for gr in gear_ratio_list:
    # Output torque = T_motor_max * gear ratio
    T_out = T_motor_max * gr
    eq_omega = find_equilibrium_omega_airfoil(T_out, omega_max=3000.0, tol=0.001)
    # Now check motor no-load limit:
    #   rod_omega <= (motor_no_load_omega / gear_ratio)
    # If eq_omega is bigger than that, it's not feasible to actually spin that fast.
    max_rod_omega_if_no_load = motor_no_load_omega / gr
    
    if eq_omega > max_rod_omega_if_no_load:
        # Not feasible: The required rod speed for equilibrium is beyond motor no-load limit
        equilibrium_rpm_list.append(float('nan'))
        feasible_list.append(False)
    else:
        feasible_list.append(True)
        # Convert eq_omega -> RPM
        eq_rpm = (eq_omega / (2.0*math.pi))*60.0
        equilibrium_rpm_list.append(eq_rpm)

print("FINAL EQUILIBRIUM SPEED CHECK (VS. NO-LOAD FEASIBILITY):")
print("GearRatio |  EquilRPM | Feasible? | Comments")
for gr, eq_rpm, feas in zip(gear_ratio_list, equilibrium_rpm_list, feasible_list):
    if feas:
        print(f"{gr:8.1f} | {eq_rpm:9.2f} |    YES    | ")
    else:
        print(f"{gr:8.1f} |       --- |    NO     | Exceeds motor no-load speed")

print("--------------------------------------------------")

# ------------------------
# 6) Plot the final RPM vs. gear ratio
# ------------------------
plt.figure(figsize=(9,5))
valid_ratios = []
valid_rpms = []
for (gr, rpm, feas) in zip(gear_ratio_list, equilibrium_rpm_list, feasible_list):
    if feas and not math.isnan(rpm):
        valid_ratios.append(gr)
        valid_rpms.append(rpm)

plt.plot(valid_ratios, valid_rpms, 'bo-', label='Feasible Equilibrium RPM')
# Mark the "not feasible" points
for (gr, rpm, feas) in zip(gear_ratio_list, equilibrium_rpm_list, feasible_list):
    if not feas:
        # Plot a red X
        plt.plot(gr, 0, 'rx', label='_nolegend_')  # put them on the baseline for a visual
        # or you could do: plt.plot(gr, 0, 'rx')
plt.title("Final Equilibrium RPM vs. Gear Ratio\n(Airfoil Drag, Checking No-load Feasibility)")
plt.xlabel("Gear Ratio (motor : rod)")
plt.ylabel("Rod RPM at Equilibrium (if feasible)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# ------------------------
# 7) For reference: plot the drag torque vs. omega
# ------------------------
omega_plot_range = range(0, 201, 10)
drag_vals = [rod_drag_torque_airfoil(o) for o in omega_plot_range]

plt.figure(figsize=(9,5))
plt.plot(omega_plot_range, drag_vals, 'r-o')
plt.title("Rod Drag Torque vs. Angular Velocity\n(Airfoil Model)")
plt.xlabel("Omega (rad/s)")
plt.ylabel("Drag Torque (N*m)")
plt.grid(True)
plt.tight_layout()
plt.show()

print("DONE. The table above shows if each gear ratio is feasible (i.e. does not exceed motor no-load speed).")
