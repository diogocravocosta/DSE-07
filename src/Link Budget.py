import numpy as np

import numpy as np

def compute_link_budget(EIRP_dBW, distance_km, frequency_GHz, G_over_T_dB, bitrate_bps, L_atm_dB=0):
    """
    EIRP_dBW: Effective isotropic radiated power [dBW]
    distance_km: Slant range to receiver [km]
    frequency_GHz: Carrier frequency [GHz]
    G_over_T_dB: Receiver G/T [dB/K]
    bitrate_bps: Bitrate [bps]
    L_atm_dB: Atmospheric loss [dB] (default = 0)
    """
    # Free Space Path Loss (FSPL)
    FSPL = 20 * np.log10(distance_km) + 20 * np.log10(frequency_GHz) + 92.45

    # Convert bitrate to dBHz
    Rb_dBHz = 10 * np.log10(bitrate_bps)

    # Compute Eb/No with atmospheric attenuation included
    Eb_No_dB = EIRP_dBW - (FSPL + L_atm_dB) + G_over_T_dB + 228.6 - Rb_dBHz

    return Eb_No_dB

# Example input values
EIRP_dBW = 50 + 0            # Transmitter EIRP [dBW]
distance_km = 600        # Distance to vehicle [km]
frequency_GHz = 2.2      # S-band
G_over_T_dB_downlink = 15         # Receiver G/T [dB/K]
G_over_T_dB_uplink = 5         # Receiver G/T [dB/K]
downlink_bitrate_bps = 10e6        # Bitrate = 1 Mbps
uplink_bitrate_bps = 0.1e6        # Bitrate = 1 Mbps
L_atm_dB = 3           # Estimated atmospheric attenuation for S-band [dB]
Eb_No_required = 3.5 #dB
# Compute link margin
Eb_No_uplink = compute_link_budget(EIRP_dBW, distance_km, frequency_GHz, G_over_T_dB_uplink, uplink_bitrate_bps, L_atm_dB)
print(f"Eb/No: {Eb_No_uplink:.2f} dB")

Eb_No_downlink = compute_link_budget(EIRP_dBW, distance_km, frequency_GHz, G_over_T_dB_downlink, downlink_bitrate_bps, L_atm_dB)
print(f"Eb/No: {Eb_No_downlink:.2f} dB")

print("Closing uplink budget: " + str(Eb_No_required - Eb_No_uplink))
print("Closing downlink budget: " + str(Eb_No_required - Eb_No_downlink))