from h2ermes_tools.propulsion.expansion_ratio_optimizer import optimize_expansion_ratio
from h2ermes_tools.propulsion.Chamber_Sizer import obtain_exit_pressure, obtain_cstar, throat_geometry, exit_geometry, chamber_geometry
from h2ermes_tools.propulsion.Mixture_Ratio_Optimizer import CEA_Obj
import numpy as np
import math
import pytest
from unittest.mock import patch
from unittest.mock import MagicMock

# Constants
g0 = 9.80665
expansion_ratios = np.array([40, 50, 60, 70, 80, 85, 90, 95, 100, 110, 120, 130, 140])
vacuum_Isp = np.array([435.5235, 439.7381, 443.0285, 445.669, 447.8579, 448.8227, 449.7163, 450.5474, 451.3234, 452.7236, 453.9782, 455.0922, 456.0993])
sea_level_Isp = np.array([274.4114, 238.3587, 201.3354, 163.7011, 125.6168, 106.4455, 87.2033, 67.8989, 48.5397, 9.6704, -29.7567, -68.9614, -108.2737])
structural_mass = 21
delta_V_vacuum = 7264.29
delta_V_sea_level = 250
payload_mass = 15000

def test_output_lengths():
    with patch("matplotlib.pyplot.show"):  # Prevent plot rendering
        filtered_ratios, total_mass_list = optimize_expansion_ratio(
            expansion_ratios, vacuum_Isp, sea_level_Isp, g0, structural_mass,
            delta_V_vacuum, delta_V_sea_level, payload_mass
        )
    
    # Check the lengths match and negative sea-level Isp values were filtered
    assert len(filtered_ratios) == len(total_mass_list)
    assert all(sea_level_Isp[sea_level_Isp >= 0])  # Only non-negative used
    assert len(filtered_ratios) == sum(sea_level_Isp >= 0)

def test_known_result():
    # Test a known input with manually calculated values
    single_exp = np.array([40])
    single_vac_isp = np.array([435.5235])
    single_sea_isp = np.array([274.4114])
    
    with patch("matplotlib.pyplot.show"):
        filtered_ratios, total_masses = optimize_expansion_ratio(
            single_exp, single_vac_isp, single_sea_isp, g0,
            structural_mass, delta_V_vacuum, delta_V_sea_level, payload_mass
        )
    
    # Recalculate manually
    mass_ratio_vac = math.exp(delta_V_vacuum / (435.5235 * g0))
    prop_mass_vac = mass_ratio_vac * (payload_mass + structural_mass) - structural_mass - payload_mass
    mass_ratio_land = math.exp(delta_V_sea_level / (274.4114 * g0))
    prop_mass_land = structural_mass * (mass_ratio_land - 1)
    expected_total = prop_mass_vac + prop_mass_land

    assert len(filtered_ratios) == 1
    assert np.isclose(total_masses[0], expected_total, rtol=1e-5)

def test_filtering_behavior():
    # Test that negative sea level ISPs are removed
    test_sea_level_Isp = np.array([100, -10, 200])
    test_vacuum_Isp = np.array([400, 410, 420])
    test_exp_ratios = np.array([10, 20, 30])
    
    with patch("matplotlib.pyplot.show"):
        filtered_ratios, _ = optimize_expansion_ratio(
            test_exp_ratios, test_vacuum_Isp, test_sea_level_Isp, g0,
            structural_mass, delta_V_vacuum, delta_V_sea_level, payload_mass
        )
    
    assert np.all(filtered_ratios == np.array([10, 30]))  # Should remove index 1

@pytest.fixture
def mock_cea():
    """Mock the CEA_Obj class and its methods."""
    with patch('h2ermes_tools.propulsion.Mixture_Ratio_Optimizer.CEA_Obj') as MockCEA:
        mock = MagicMock()
        mock.get_Cstar.return_value = 1580
        mock.get_PcOvPe.return_value = 120
        mock.get_Chamber_MolWt_gamma.return_value = [0, 0, 1.22]  # gamma = 1.22
        MockCEA.return_value = mock
        yield mock


def test_obtain_cstar(mock_cea):
    Pc = 6e6
    MR = 6.0
    cstar = obtain_cstar(Pc, MR)
    assert cstar == 1580
    mock_cea.get_Cstar.assert_called_once_with(Pc, MR)


def test_obtain_exit_pressure(mock_cea):
    Pc = 6e6
    MR = 6.0
    eps = 80
    p_exit = obtain_exit_pressure(Pc, MR, eps)
    expected = Pc / 120
    assert np.isclose(p_exit, expected)


def test_throat_geometry(mock_cea):
    thrust = 66700
    Isp = 450
    Pc = 6e6
    MR = 6.0

    area, diameter = throat_geometry(thrust, Isp, Pc, MR)
    mass_flow = thrust / (9.81 * Isp)
    expected_area = mass_flow * 1580 / Pc
    expected_diameter = np.sqrt(4 * expected_area / np.pi)
    
    assert np.isclose(area, expected_area)
    assert np.isclose(diameter, expected_diameter)


def test_exit_geometry(mock_cea):
    thrust = 66700
    Isp = 450
    Pc = 6e6
    MR = 6.0
    eps = 80

    area, diameter = exit_geometry(thrust, Isp, Pc, MR, eps)
    throat_area = throat_geometry(thrust, Isp, Pc, MR)[0]
    expected_area = throat_area * eps
    expected_diameter = np.sqrt(4 * expected_area / np.pi)

    assert np.isclose(area, expected_area)
    assert np.isclose(diameter, expected_diameter)


def test_chamber_geometry(mock_cea):
    thrust = 66700
    Isp = 450
    Pc = 6e6
    MR = 6.0
    eps = 80
    L_star = 0.76

    values = chamber_geometry(thrust, Isp, Pc, MR, eps, L_star)
    combustion_length, D_c, a_t, D_t, D_e, contraction_ratio, area = values

    assert all(isinstance(v, float) for v in values)
    assert contraction_ratio > 1
    assert D_c > 0
    assert a_t > 0
    assert area > a_t



from h2ermes_tools.propulsion.feed_system_sizing import obtain_size_propellant_channel
def test_typical_case():
    mass_flow = 15.2  # kg/s
    density = 70.85   # kg/m^3
    velocity = 10     # m/s
    k = 1.2           # head loss coefficient

    area, diameter = obtain_size_propellant_channel(mass_flow, density, velocity, k)

    expected_area = mass_flow / (density * velocity)
    expected_diameter = 2 * np.sqrt(expected_area / np.pi)

    assert np.isclose(area, expected_area, rtol=1e-6)
    assert np.isclose(diameter, expected_diameter, rtol=1e-6)


def test_zero_mass_flow():
    area, diameter = obtain_size_propellant_channel(0, 70.85, 10, 1.2)
    assert area == 0
    assert diameter == 0


def test_high_velocity():
    area, diameter = obtain_size_propellant_channel(15.2, 70.85, 100, 1.2)
    assert area < 0.01
    assert diameter < 0.2


def test_negative_inputs():
    # These don't raise errors in your code, but you might want to handle them in production
    area, diameter = obtain_size_propellant_channel(-15.2, 70.85, 10, 1.2)
    assert area < 0
    assert diameter > 0  # square root of positive area still gives valid diameter (complex cases are avoided here)

    # Negative density
    area, diameter = obtain_size_propellant_channel(15.2, -70.85, 10, 1.2)
    assert area < 0

    # Negative velocity
    area, diameter = obtain_size_propellant_channel(15.2, 70.85, -10, 1.2)
    assert area < 0


def test_extremely_small_mass_flow():
    area, diameter = obtain_size_propellant_channel(1e-6, 70.85, 10, 1.2)
    assert area > 0
    assert diameter > 0
    assert diameter < 0.01


from h2ermes_tools.propulsion.Mixture_Ratio_Optimizer import eps_for_isp, calculate_propmass, Mixture_Ratio_Optimizer

C = CEA_Obj(
    oxName='LOX', fuelName='LH2',
    pressure_units='Pa', cstar_units='m/s',
    temperature_units='K', sonic_velocity_units='m/s',
    enthalpy_units='J/kg', density_units='kg/m^3',
    specific_heat_units='J/kg-K', viscosity_units='millipoise',
    thermal_cond_units='mcal/cm-K-s'
)

def test_eps_for_isp_typical():
    eps, isp, correction = eps_for_isp(MR=6, pc=6.1e6, Isp_desired=450)
    assert 50 < eps < 200, "Expansion ratio out of typical range"
    assert 440 < isp < 470, "Corrected Isp should be near target"
    assert 0 < correction < 10, "Correction factor % should be reasonable"

def test_eps_for_isp_high_Isp():
    eps, isp, correction = eps_for_isp(MR=6, pc=6.1e6, Isp_desired=470)
    assert eps > 100, "Higher Isp requires higher expansion ratio"

def test_calculate_propmass_basic():
    result = calculate_propmass(struct_ratio=0.1, Isp=450, deltaV=9000, payload=1000)
    assert result > 0
    assert isinstance(result, float)

def test_calculate_propmass_edge_zero_payload():
    result = calculate_propmass(struct_ratio=0.1, Isp=450, deltaV=9000, payload=0)
    assert result == 0

def test_Mixture_Ratio_Optimizer_typical():
    mr, isp, mass = Mixture_Ratio_Optimizer(
        pc=6.1e6, struct_ratio=0.12, deltaV=7500, payload=10000,
        OF_MIN=4, OF_MAX=8, eps=80
    )
    assert 4 <= mr <= 8
    assert 400 <= isp <= 480
    assert mass > 0

def test_Mixture_Ratio_Optimizer_increasing_payload():
    _, _, mass1 = Mixture_Ratio_Optimizer(payload=10000)
    _, _, mass2 = Mixture_Ratio_Optimizer(payload=20000)
    assert mass2 > mass1, "Higher payload should require more propellant"



