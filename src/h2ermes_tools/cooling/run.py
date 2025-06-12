from pyfluids import Fluid, FluidsList, Input
from h2ermes_tools.cooling.coolant import Coolant
from h2ermes_tools.cooling.channel import RectangularChannel
from h2ermes_tools.variables import coolant_inlet_pressure, coolant_inlet_temperature
from h2ermes_tools.cooling.heat_transfer import HeatShield
from h2ermes_tools.cooling.material import SS310
import numpy as np

if __name__ == "__main__":
    # --- Channel and Coolant Setup ---
    channel_length = 5.0  # [m]
    wall_thickness = 5e-3  # [m]
    num_nodes = 3
    # Initial coolant state
    hydrogen = Fluid(FluidsList.Hydrogen).with_state(
        Input.pressure(coolant_inlet_pressure.value),
        Input.temperature(coolant_inlet_temperature.value),
    )
    mass_flow = 0.0001
    coolant = Coolant(
        fluid=hydrogen,
        channel=RectangularChannel(
            length=channel_length,
            width=5e-3,  # [m]
            height=5e-3,  # [m]
            roughness=1e-5,  # [m]
        ),
        mass_flow=mass_flow,
    )

    # --- Simulation Control ---
    do_animation = False
    incident_heat_flux_profile = np.array(
        [
            [0, 100_000],
            [150, 666_000],
            [450, 666_000],
            [900, 0],
        ]
    )
    initial_temperature = 300.0
    total_time = 900.0
    hs = HeatShield(
        wall_thickness=wall_thickness,
        num_nodes=num_nodes,
        incident_heat_flux=incident_heat_flux_profile,
        initial_temperature=initial_temperature,
        total_time=total_time,
        coolant=coolant,
        material=SS310,
        heat_shield_diameter=10,
        sphere_height= 2.0,
    )

    print(hs.estimate_heat_shield_mass(h=2))
    print(hs.coolant.channel.width)
    print(hs.coolant.channel.height)
    print(hs.wall_thickness)

    # all_temperatures, all_times, x, fourier_number, time_step, node_spacing = (
    #     hs.run_simulation()
    # )

    # if do_animation:
    #     HeatShield.animate_profiles(
    #         all_temperatures,
    #         all_times,
    #         x,
    #         interval=2,
    #     )
    # else:
    #     HeatShield.plot_profiles(
    #         all_temperatures,
    #         all_times,
    #         x,
    #         num_profiles=10,
    #     )

    # # --- Space-Time Simulation (Corrected Marching) ---
    # num_segments = 4
    # segment_length = channel_length / (num_segments + 1)
    # n_timesteps = 10000  # manually set number of time steps
    # wall_states = [np.full(num_nodes, initial_temperature, dtype=float) for _ in range(num_segments)]
    # # Initialize coolant state for each segment (first segment gets inlet, others get copies)
    # coolant_states = []
    # for _ in range(num_segments):
    #     fluid_new = Fluid(FluidsList.Hydrogen).with_state(
    #         Input.pressure(coolant.fluid.pressure),
    #         Input.temperature(coolant.fluid.temperature),
    #     )
    #     coolant_states.append(Coolant(
    #         fluid=fluid_new,
    #         channel=coolant.channel,
    #         mass_flow=coolant.mass_flow,
    #     ))
    # all_wall_histories = [[] for _ in range(num_segments)]
    # all_coolant_histories = [[] for _ in range(num_segments)]
    # hs_temp = HeatShield(
    #     wall_thickness=wall_thickness,
    #     num_nodes=num_nodes,
    #     incident_heat_flux=incident_heat_flux_profile,
    #     initial_temperature=initial_temperature,
    #     total_time=total_time,
    #     coolant=coolant,
    #     material=SS310,
    #     heat_shield_diameter=10,
    # )
    # density = hs_temp.heat_shield_density
    # specific_heat = hs_temp.specific_heat
    # thermal_conductivity = hs_temp.thermal_conductivity
    # T_ref = initial_temperature
    # cp_ref = specific_heat(T_ref) if callable(specific_heat) else specific_heat
    # k_ref = thermal_conductivity(T_ref) if callable(thermal_conductivity) else thermal_conductivity
    # thermal_diffusivity = k_ref / (density * cp_ref)
    # node_spacing = wall_thickness / (num_nodes - 1)
    # max_dt = 0.5 * node_spacing**2 / thermal_diffusivity
    # dt = 0.05 * max_dt  # manually adjusted to ensure stability
    # times = [0.0]
    # for t_idx in range(n_timesteps):
    #     t_now = t_idx * dt
    #     times.append(t_now)
    #     # At the start of each time step, propagate coolant states from previous time step
    #     coolant_in_segments = []
    #     for seg in range(num_segments):
    #         if t_idx == 0:
    #             # First time step: all segments get inlet coolant
    #             coolant_in_segments.append(coolant_states[seg])
    #         elif seg == 0:
    #             # First segment: propagate from previous time step of same segment
    #             prev_temp, prev_pres = all_coolant_histories[seg][-1]
    #             fluid_new = Fluid(FluidsList.Hydrogen).with_state(
    #                 Input.pressure(prev_pres),
    #                 Input.temperature(prev_temp),
    #             )
    #             coolant_in_segments.append(Coolant(
    #                 fluid=fluid_new,
    #                 channel=coolant.channel,
    #                 mass_flow=coolant.mass_flow,
    #             ))
    #         else:
    #             # Other segments: get coolant state from previous segment (last time step)
    #             prev_temp, prev_pres = all_coolant_histories[seg-1][-1]
    #             fluid_new = Fluid(FluidsList.Hydrogen).with_state(
    #                 Input.pressure(prev_pres),
    #                 Input.temperature(prev_temp),
    #             )
    #             coolant_in_segments.append(Coolant(
    #                 fluid=fluid_new,
    #                 channel=coolant.channel,
    #                 mass_flow=coolant.mass_flow,
    #             ))
    #     for seg in range(num_segments):
    #         hs_seg = HeatShield(
    #             wall_thickness=wall_thickness,
    #             num_nodes=num_nodes,
    #             incident_heat_flux=incident_heat_flux_profile,
    #             initial_temperature=wall_states[seg][0],
    #             total_time=dt,
    #             coolant=coolant_in_segments[seg],
    #             material=SS310,
    #             heat_shield_diameter=10,
    #         )
    #         hs_seg_wall = wall_states[seg].copy()
    #         density = hs_seg.heat_shield_density
    #         specific_heat = hs_seg.specific_heat
    #         thermal_conductivity = hs_seg.thermal_conductivity
    #         stefan_boltzmann = hs_seg.stefan_boltzmann
    #         node_spacing = wall_thickness / (num_nodes - 1)
    #         temperature = hs_seg_wall
    #         temperature_next = temperature.copy()
    #         times_hf = incident_heat_flux_profile[:, 0]
    #         values_hf = incident_heat_flux_profile[:, 1]
    #         q_incident = np.interp(t_now, times_hf, values_hf)
    #         cp_nodes = specific_heat(temperature)
    #         k_nodes = thermal_conductivity(temperature)
    #         for idx in range(1, num_nodes - 1):
    #             k_left = k_nodes[idx - 1]
    #             k_right = k_nodes[idx]
    #             k_avg = 0.5 * (k_left + k_right)
    #             temperature_next[idx] = temperature[idx] + (
    #                 dt * k_avg / (density * cp_nodes[idx] * node_spacing**2) * (temperature[idx + 1] - 2 * temperature[idx] + temperature[idx - 1])
    #             )
    #         radiative_loss = SS310.emissivity * stefan_boltzmann * temperature[0] ** 4
    #         k0 = k_nodes[0]
    #         cp0 = cp_nodes[0]
    #         net_heat_flux_0 = (
    #             q_incident - radiative_loss - k0 * (temperature[0] - temperature[1]) / node_spacing
    #         )
    #         control_volume_0 = node_spacing / 2.0
    #         temperature_next[0] = temperature[0] + dt * net_heat_flux_0 / (
    #             density * cp0 * control_volume_0
    #         )
    #         kN = k_nodes[-1]
    #         cpN = cp_nodes[-1]
    #         cold_wall_temp = temperature[num_nodes - 1]
    #         T_bulk = coolant_in_segments[seg].fluid.temperature
    #         h_cool = coolant_in_segments[seg].get_heat_transfer_coefficient()
    #         area = coolant_in_segments[seg].channel.get_contact_area(segment_length)
    #         cooling_flux = h_cool * (cold_wall_temp - T_bulk)
    #         net_heat_flux_last = (
    #             kN * (temperature[num_nodes - 2] - temperature[num_nodes - 1]) / node_spacing - cooling_flux
    #         )
    #         control_volume_last = node_spacing / 2.0
    #         temperature_next[num_nodes - 1] = temperature[
    #             num_nodes - 1
    #         ] + dt * net_heat_flux_last / (density * cpN * control_volume_last)
    #         energy_to_coolant = cooling_flux * area * dt
    #         coolant_in_segments[seg].add_energy(energy_to_coolant)
    #         dP = coolant_in_segments[seg].calculate_pressure_drop(segment_length)
    #         coolant_in_segments[seg].reduce_pressure(dP)
    #         wall_states[seg] = temperature_next.copy()
    #         all_wall_histories[seg].append(temperature_next.copy())
    #         all_coolant_histories[seg].append(
    #             (coolant_in_segments[seg].fluid.temperature, coolant_in_segments[seg].fluid.pressure)
    #         )
    # # Print coolant state at last time step for each segment
    # for i in range(num_segments):
    #     temp, pres = all_coolant_histories[i][-1]
    #     print(
    #         f"Segment {i + 1}: Coolant Pressure = {pres * 1e-5:.2f} bar, Coolant Temperature = {temp:.2f} K"
    #     )
    # wall_hist = np.array(all_wall_histories[-1])
    # times_arr = np.arange(n_timesteps) * dt
    # if do_animation:
    #     HeatShield.animate_profiles(
    #         wall_hist,
    #         times_arr,
    #         x,
    #         interval=2,
    #     )
    # else:
    #     HeatShield.plot_profiles(
    #         wall_hist,
    #         times_arr,
    #         x,
    #         num_profiles=10,
    #     )

    # # TODO: doesn't work right now
    # # hs.print_summary(all_temperatures, x[1] - x[0], time_step, fourier_number)
