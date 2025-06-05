from h2ermes_tools.n2_chart_generator import plot_n2

if __name__ == "__main__":
    elements_io = {
        "Cool heat shield": {
            "inputs": [],
            "outputs": [("Pump coolant", "α"), ("Sustain re-entry\n loads", "β")]
        },
        "Pump coolant": {
            "inputs": [("Cool heat shield", "α")],
            "outputs": [("Sustain re-entry\n loads", "γ")]
        },
        "Sustain re-entry\n loads": {
            "inputs": [("Cool heat shield", "β"), ("Pump coolant", "γ")],
            "outputs": []
        },
        "Provide power": {
            "inputs": [],
            "outputs": [("Cool heat shield", "δ"), ("Pump coolant", "ε")]
        },
        "Provide structural\n support": {
            "inputs": [("Cool heat shield", "ζ")],
            "outputs": [("Sustain re-entry\n loads", "η")]
        },
        "Provide rigidity": {
            "inputs": [("Cool heat shield", "θ"), ("Pump coolant", "ι")],
            "outputs": []
        },
        "Provide cooling": {
            "inputs": [],
            "outputs": [("Cool heat shield", "κ"), ("Pump coolant", "λ")]
        },
        "Provide thermal\n insulation": {
            "inputs": [("Cool heat shield", "μ")],
            "outputs": [("Sustain re-entry\n loads", "ν")]
        },
        "Provide thermal\n protection": {
            "inputs": [("Cool heat shield", "ξ"), ("Pump coolant", "ο")],
            "outputs": []
        },
        "Provide\n aerodynamics": {
            "inputs": [],
            "outputs": [("Cool heat shield", "π"), ("Pump coolant", "ρ")]
        },
    }

    # Plot the N2 chart
    plot_n2(elements_io, title="N2 Chart Example", save=True)
