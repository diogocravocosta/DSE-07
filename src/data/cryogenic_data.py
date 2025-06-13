import math
import numpy as np
import matplotlib.pyplot as plt

# --- Data Structure for Material Properties ---
MATERIALS = {
    "SS304": {
        "Thermal Conductivity": [
            {
                "coefficients": [
                    -1.4087,
                    1.3982,
                    0.2543,
                    -0.6260,
                    0.2334,
                    0.4256,
                    -0.4658,
                    0.1650,
                    -0.0199,
                ],
                "range": (1, 300),
                "unit": "W/(m·K)",
            }
        ],
        "Specific Heat": [
            {
                "coefficients": [
                    22.0061,
                    -127.5528,
                    303.647,
                    -381.0098,
                    274.0328,
                    -112.9212,
                    24.7593,
                    -2.239153,
                    0,
                ],
                "range": (4, 300),
                "unit": "J/(kg·K)",
            }
        ],
    },
    "SS304L": {
        "Thermal Conductivity": [
            {
                "coefficients": [
                    -1.4087,
                    1.3982,
                    0.2543,
                    -0.6260,
                    0.2334,
                    0.4256,
                    -0.4658,
                    0.1650,
                    -0.0199,
                ],
                "range": (1, 300),
                "unit": "W/(m·K)",
            }
        ],
        "Specific Heat": [
            {
                "coefficients": [
                    -351.51,
                    3123.695,
                    -12017.28,
                    26143.99,
                    -35176.33,
                    29981.75,
                    -15812.78,
                    4719.64,
                    -610.515,
                ],
                "range": (4, 23),
                "unit": "J/(kg·K)",
            }
        ],
    },
    "SS310": {
        "Thermal Conductivity": [
            {
                "coefficients": [
                    -0.81907,
                    -2.1967,
                    9.1059,
                    -13.078,
                    10.853,
                    -5.1269,
                    1.2583,
                    -0.12395,
                    0,
                ],
                "range": (1, 300),
                "unit": "W/(m·K)",
            }
        ],
        "Specific Heat": [
            {  # Range 1 for Specific Heat
                "coefficients": [
                    20.694,
                    -171.007,
                    600.6256,
                    -1162.748,
                    1361.931,
                    -986.2934,
                    430.093,
                    -102.85,
                    10.275,
                ],
                "range": (4, 47),
                "unit": "J/(kg·K)",
            },
            {  # Range 2 for Specific Heat
                "coefficients": [
                    -2755.63,
                    9704.831,
                    -14618.36,
                    12202.74,
                    -6092.339,
                    1818.555,
                    -300.458,
                    21.1942,
                    0,
                ],
                "range": (47, 300),
                "unit": "J/(kg·K)",
            },
        ],
    },
}


def _calculate_polynomial(T, coefficients):
    """Internal function to solve the log-polynomial equation."""
    log10_T = math.log10(T)
    log10_y = sum(c * (log10_T**i) for i, c in enumerate(coefficients))
    return 10**log10_y


def get_property(material_name, property_name, T):
    """
    Calculates a material property for a given temperature.

    Args:
        material_name (str): The name of the material (e.g., 'SS310').
        property_name (str): The name of the property (e.g., 'Thermal Conductivity').
        T (float): The temperature in Kelvin.

    Returns:
        float: The calculated material property value.

    Raises:
        ValueError: If the material, property, or temperature is not valid.
    """
    if material_name not in MATERIALS or property_name not in MATERIALS[material_name]:
        raise ValueError(
            f"Property '{property_name}' for material '{material_name}' not found."
        )

    property_data_list = MATERIALS[material_name][property_name]

    for prop_data in property_data_list:
        min_temp, max_temp = prop_data["range"]
        if min_temp <= T <= max_temp:
            return _calculate_polynomial(T, prop_data["coefficients"])

    # If no valid range was found after checking all lists
    valid_ranges = " or ".join(
        [f"{p['range'][0]}-{p['range'][1]}" for p in property_data_list]
    )
    raise ValueError(
        f"Temperature {T} K is outside the valid range(s) {valid_ranges} K for {property_name}."
    )


def plot_properties(properties_to_plot, title_override=None, num_points=500):
    """
    Generates and displays a plot for one or more material properties.

    Args:
        properties_to_plot (list): A list of tuples, where each tuple contains
                                   (material_name, property_name).
        title_override (str, optional): A custom title for the plot.
        num_points (int): The number of data points to generate for the plot.
    """
    fig, ax = plt.subplots()

    y_unit = ""
    plot_title_parts = []

    for material_name, property_name in properties_to_plot:
        property_data_list = MATERIALS[material_name][property_name]

        # Determine a common unit for the y-axis
        if not y_unit:
            y_unit = property_data_list[0]["unit"]
        elif y_unit != property_data_list[0]["unit"]:
            print(
                f"Warning: Plotting properties with different units ({y_unit} and {property_data_list[0]['unit']})."
            )

        if f"{material_name}" not in plot_title_parts:
            plot_title_parts.append(f"{material_name}")

        for prop_data in property_data_list:
            min_temp, max_temp = prop_data["range"]
            temp_range = np.linspace(min_temp, max_temp, num_points)
            values = [
                _calculate_polynomial(T, prop_data["coefficients"]) for T in temp_range
            ]

            label = f"{material_name}"
            if len(property_data_list) > 1:
                label += f" ({min_temp}-{max_temp} K)"

            ax.plot(temp_range, values, label=label)

    # Determine plot title
    if title_override:
        final_title = title_override
    else:
        # Use the property name as the title if all items are the same property
        main_prop_name = properties_to_plot[0][1]
        if all(p[1] == main_prop_name for p in properties_to_plot):
            final_title = f"{main_prop_name} of {', '.join(plot_title_parts)}"
        else:
            final_title = f"Cryogenic Properties of {', '.join(plot_title_parts)}"

    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(f"Value ({y_unit})")
    ax.set_title(final_title)
    ax.legend()
    ax.grid(True)
    ax.set_yscale(
        "log"
    )  # Use a log scale for better visualization of wide-ranging values
    plt.show()


if __name__ == "__main__":
    # Example 1: Calculate properties for the new SS304 at T = 100 K
    temp = 100.0
    print(f"Calculating properties for SS304 at T = {temp} K...")
    try:
        tc = get_property("SS304", "Thermal Conductivity", temp)
        print(f"- Thermal Conductivity: {tc:,.4f} W/(m·K)")

        sh = get_property("SS304", "Specific Heat", temp)
        print(f"- Specific Heat: {sh:,.4f} J/(kg·K)")
    except ValueError as e:
        print(e)

    print("Generating Comparative Plots...")

    # Example 2: Plot a comparison of Thermal Conductivity for all steels
    print("Plotting Thermal Conductivity for all steels...")
    plot_properties(
        [
            ("SS304", "Thermal Conductivity"),
            ("SS304L", "Thermal Conductivity"),
            ("SS310", "Thermal Conductivity"),
        ]
    )

    # Example 3: Plot a comparison of Specific Heat for all steels
    print("Plotting Specific Heat for all steels...")
    plot_properties(
        [
            ("SS304", "Specific Heat"),
            ("SS304L", "Specific Heat"),  # Note: this has a very narrow range
            ("SS310", "Specific Heat"),
        ]
    )

    # Note: specific heat data is readily available from Granta
    # This script is only useful for thermal conductivity at cryogenic temperatures
