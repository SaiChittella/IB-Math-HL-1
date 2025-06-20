import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from tqdm import tqdm # For progress bar
import calendar # To get number of days in a month more accurately

# Constants
deg2rad = np.pi / 180
rad2deg = 180 / np.pi

# --- 1. Declination angle (degrees)
def declination_angle(n):
    return 23.45 * np.sin(deg2rad * (360 / 365 * (284 + n)))

# --- 2. Hour angle at time t (solar time, in hours)
def hour_angle(solar_time):
    return 15 * (solar_time - 12)  # degrees

# --- 3. Solar altitude angle (h in degrees)
def solar_altitude(phi, delta, omega):
    phi, delta, omega = map(np.radians, [phi, delta, omega])
    sin_h = np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.cos(omega)
    return np.degrees(np.arcsin(sin_h))

# --- 4. Solar azimuth angle (degrees)
def solar_azimuth(phi, delta, omega, h):
    # all in degrees
    phi, delta, omega, h = map(np.radians, [phi, delta, omega, h])
    
    cos_A = (np.sin(delta) * np.cos(phi) - np.cos(delta) * np.sin(phi) * np.cos(omega)) / np.cos(h)
    cos_A = np.clip(cos_A, -1, 1)
    A = np.arccos(cos_A)
    
    if omega > 0:
        A = 2 * np.pi - A  # afternoon
    return np.degrees(A)

# --- 5. Sunset hour angle (degrees)
def sunset_hour_angle(phi, delta):
    phi, delta = map(np.radians, [phi, delta])
    cos_omega_ss = -np.tan(phi) * np.tan(delta)
    cos_omega_ss = np.clip(cos_omega_ss, -1, 1)
    return np.degrees(np.arccos(cos_omega_ss))

# --- 6. Sunrise and Sunset solar times
def sunrise_sunset_times(phi, delta):
    omega_ss = sunset_hour_angle(phi, delta)
    daylight_hours = 2 * omega_ss / 15
    sunrise = 12 - omega_ss / 15
    sunset = 12 + omega_ss / 15
    return sunrise, sunset, daylight_hours

# --- 7. Cosine of incidence angle (between sun and panel)
def cos_theta_i(h, A_sun, beta, A_panel):
    h = np.radians(h)
    A_sun = np.radians(A_sun)
    beta = np.radians(beta)
    A_panel = np.radians(A_panel)
    return (np.sin(h) * np.cos(beta) +
            np.cos(h) * np.sin(beta) * np.cos(A_sun - A_panel))

# --- 8. Power from beam at time t (W/m^2), assuming constant I_b
def power_beam(t, phi, delta, beta, A_panel, I_b):
    omega = hour_angle(t)
    h = solar_altitude(phi, delta, omega)
    if h <= 0:
        return 0
    A_sun = solar_azimuth(phi, delta, omega, h)
    
    cos_theta = cos_theta_i(h, A_sun, beta, A_panel)
    power = I_b * max(0, cos_theta) # Ensure power is not negative
    return power

# --- 9. Total beam energy over a day (Wh/m^2)
def total_daily_beam_energy(phi, n, beta, A_panel, I_b=1000):
    delta = declination_angle(n)
    sunrise, sunset, _ = sunrise_sunset_times(phi, delta)
    
    # If the sun never rises (e.g., polar winter)
    if sunrise >= sunset:
        return 0 
        
    integrand = lambda t: power_beam(t, phi, delta, beta, A_panel, I_b)
    
    # Using a robust integration approach for potential discontinuities
    energy, _ = quad(integrand, sunrise, sunset, epsabs=1e-2, limit=100) 
    
    return energy

# --- Function: Total Energy over a period ---
def total_energy_over_period(phi, start_day, end_day, beta, A_panel, I_b=1000):
    total_energy = 0
    for n_day in range(start_day, end_day + 1):
        total_energy += total_daily_beam_energy(phi, n_day, beta, A_panel, I_b)
    return total_energy

# --- Function: Get days for a given month ---
def get_month_days(month_index, year=2024): # Added year for calendar.monthrange
    # month_index: 1 for Jan, 2 for Feb, ..., 12 for Dec
    month_names = ["", "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
    
    num_days_in_month = calendar.monthrange(year, month_index)[1]
    
    # Calculate start day of the year for the given month
    start_day_of_year = 0
    for m in range(1, month_index):
        start_day_of_year += calendar.monthrange(year, m)[1]
    
    start_day = start_day_of_year + 1
    end_day = start_day_of_year + num_days_in_month
    
    return start_day, end_day, month_names[month_index], num_days_in_month


# --- Main Parameters ---
phi = 47.633  # Sammamish latitude
A_panel = 180 # facing due south
I_b = 1000    # beam irradiance (W/m²)
year_for_calculation = 2024 # Use a non-leap year for consistency with 365, or adjust declination if using leap year

# --- Focus on July ---
selected_month_index = 7 # July
selected_month_start_day, selected_month_end_day, selected_month_name, num_days_in_july = get_month_days(selected_month_index, year_for_calculation)
selected_day_for_plot = selected_month_start_day + (num_days_in_july // 2) # Mid-month day for daily plots

print(f"--- Analysis Focused on {selected_month_name} ({num_days_in_july} days, Representative Day: {selected_day_for_plot}) ---")

# Debug values at solar noon for the selected mid-month day
delta_selected_day = declination_angle(selected_day_for_plot)
omega = 0  # solar noon
h_selected_day = solar_altitude(phi, delta_selected_day, omega)
print(f"🕛 Solar noon altitude on day {selected_day_for_plot}: {round(h_selected_day, 2)} degrees")
print(f"📅 Declination on day {selected_day_for_plot}: {round(delta_selected_day, 2)} degrees")
sunrise, sunset, daylight_hours = sunrise_sunset_times(phi, delta_selected_day)
print(f"🌅 Sunrise: {round(sunrise, 2)}, Sunset: {round(sunset, 2)}, Daylight Hours: {round(daylight_hours, 2)}")

# --- Daily Plots for the Selected July Day (Figures 1-4) ---
times = np.linspace(sunrise, sunset, 200)

altitudes = [solar_altitude(phi, delta_selected_day, hour_angle(t)) for t in times]
azimuths = [solar_azimuth(phi, delta_selected_day, hour_angle(t), solar_altitude(phi, delta_selected_day, hour_angle(t))) for t in times]

# Define specific tilts for the daily power/cumulative plot (within your new range)
betas_for_daily_plot = [5, 10, 20] # Using specific tilts from the new range
power_curves = {b: [power_beam(t, phi, delta_selected_day, b, A_panel, I_b) for t in times] for b in betas_for_daily_plot}
dt = (times[1] - times[0])
cumulative_energy_curves_daily = {b: np.cumsum(power_curves[b]) * dt for b in betas_for_daily_plot}

figs = []

# 1. Solar Altitude vs Time
fig1, ax1 = plt.subplots(figsize=(10, 4))
ax1.plot(times, altitudes)
ax1.set_title(f"Solar Altitude vs. Time (Day {selected_day_for_plot}, {selected_month_name})")
ax1.set_xlabel("Time (hours)")
ax1.set_ylabel("Altitude Angle (°)")
ax1.grid(True)
figs.append(fig1)

# 2. Solar Azimuth vs Time
fig2, ax2 = plt.subplots(figsize=(10, 4))
ax2.plot(times, azimuths)
ax2.set_title(f"Solar Azimuth vs. Time (Day {selected_day_for_plot}, {selected_month_name})")
ax2.set_xlabel("Time (hours)")
ax2.set_ylabel("Azimuth Angle (°)")
ax2.grid(True)
figs.append(fig2)

# 3. Instantaneous Power Output for Different Tilts (for the selected day)
fig3, ax3 = plt.subplots(figsize=(10, 4))
for beta, powers in power_curves.items():
    ax3.plot(times, powers, label=f"Tilt {beta}°")
ax3.set_title(f"Instantaneous Power Output vs. Time (Day {selected_day_for_plot}, {selected_month_name})")
ax3.set_xlabel("Time (hours)")
ax3.set_ylabel("Power (W/m²)")
ax3.legend()
ax3.grid(True)
figs.append(fig3)

# 4. Cumulative Energy vs. Time for Different Tilts (for the selected day) - your requested plot
fig4, ax4 = plt.subplots(figsize=(10, 4))
for beta, energies in cumulative_energy_curves_daily.items():
    ax4.plot(times, energies, label=f"Tilt {beta}°")
ax4.set_title(f"Cumulative Energy vs. Time (Day {selected_day_for_plot}, {selected_month_name})")
ax4.set_xlabel("Time (hours)")
ax4.set_ylabel("Cumulative Energy (Wh/m²)")
ax4.legend()
ax4.grid(True)
figs.append(fig4)


# --- Total Monthly Energy for July vs. Tilt Angle (Figure 5) ---
# New specified tilt range: 15 to 25 degrees with 0.5 degree interval
tilt_angles_for_monthly_optimization = np.arange(10, 12.001, 0.001) 
monthly_energies_july = []

print(f"\n--- Calculating Total Energy for {selected_month_name} at Various Tilts ({tilt_angles_for_monthly_optimization[0]}-{tilt_angles_for_monthly_optimization[-1]} deg) ---")
for tilt in tqdm(tilt_angles_for_monthly_optimization, desc=f"Calculating for {selected_month_name}"):
    monthly_e = total_energy_over_period(phi, selected_month_start_day, selected_month_end_day, tilt, A_panel, I_b)
    monthly_energies_july.append(monthly_e)

# 5. Total Energy for July vs. Tilt Angle
fig5, ax5 = plt.subplots(figsize=(10, 6))
ax5.plot(tilt_angles_for_monthly_optimization, monthly_energies_july, marker='o')
ax5.set_title(f"Total Energy for {selected_month_name} vs. Panel Tilt Angle (Phi={phi}°)")
ax5.set_xlabel("Panel Tilt Angle (°)")
ax5.set_ylabel("Total Monthly Energy (Wh/m²)")
ax5.grid(True)
ax5.set_xticks(np.arange(5, 15.1, 1)) # Show every 1 degree on x-axis
figs.append(fig5)

# Find the optimal tilt for July
optimal_tilt_index = np.argmax(monthly_energies_july)
optimal_tilt_angle_july = tilt_angles_for_monthly_optimization[optimal_tilt_index]
max_monthly_energy_july = monthly_energies_july[optimal_tilt_index]
print(f"Optimal fixed tilt for {selected_month_name}: {optimal_tilt_angle_july}° with {round(max_monthly_energy_july, 1)} Wh/m²")


# --- NEW GRAPHS ---

# 6. Daily Energy Production over the Month of July for a Fixed Tilt (Optimal July Tilt)
print(f"\n--- Calculating Daily Energy for {selected_month_name} at Optimal Tilt ({optimal_tilt_angle_july}°) ---")
daily_energies_july_optimal_tilt = []
days_of_july = range(selected_month_start_day, selected_month_end_day + 1)
for n_day in tqdm(days_of_july, desc=f"Daily energy for {selected_month_name}"):
    daily_e = total_daily_beam_energy(phi, n_day, optimal_tilt_angle_july, A_panel, I_b)
    daily_energies_july_optimal_tilt.append(daily_e)

fig6, ax6 = plt.subplots(figsize=(12, 6))
# Plot against day of month (1 to 31)
ax6.bar(range(1, num_days_in_july + 1), daily_energies_july_optimal_tilt, color='skyblue')
ax6.set_title(f"Daily Energy Production in {selected_month_name} at Optimal Tilt ({optimal_tilt_angle_july}°, Phi={phi}°)")
ax6.set_xlabel(f"Day of {selected_month_name}")
ax6.set_ylabel("Daily Energy (Wh/m²)")
ax6.set_xticks(np.arange(1, num_days_in_july + 1, 2)) # Show every other day
ax6.grid(axis='y')
figs.append(fig6)


# 7. Power Output Comparison between Optimal and a Suboptimal Tilt (for a single July day)
# We'll use the mid-month day (selected_day_for_plot) for this
# Compare optimal_tilt_angle_july with a suboptimal angle, e.g., 0 (flat) or 40 (steeper)
suboptimal_tilt_1 = 0 # flat panel
suboptimal_tilt_2 = 40 # a steeper panel

power_optimal = [power_beam(t, phi, delta_selected_day, optimal_tilt_angle_july, A_panel, I_b) for t in times]
power_suboptimal_1 = [power_beam(t, phi, delta_selected_day, suboptimal_tilt_1, A_panel, I_b) for t in times]
power_suboptimal_2 = [power_beam(t, phi, delta_selected_day, suboptimal_tilt_2, A_panel, I_b) for t in times]




plt.tight_layout()
plt.show()