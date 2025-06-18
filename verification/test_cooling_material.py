import numpy as np
from h2ermes_tools.cooling.material import Material, SS310, SS304, Ti6Al4V, SS316L

def test_material_basic_properties():
    assert isinstance(SS310.density, (int, float))
    assert isinstance(SS304.density, (int, float))
    assert isinstance(Ti6Al4V.density, (int, float))
    assert isinstance(SS316L.density, (int, float))
    assert SS310.density > 0
    assert SS304.density > 0
    assert Ti6Al4V.density > 0
    assert SS316L.density > 0

def test_material_specific_heat_callable():
    T = np.array([300, 500, 1000])
    cp = SS310.specific_heat(T)
    assert np.all(cp > 0)
    cp2 = SS304.specific_heat(T)
    assert np.all(cp2 > 0)

def test_material_thermal_conductivity_callable():
    T = np.array([300, 500, 1000])
    k = SS310.thermal_conductivity(T)
    assert np.all(k > 0)
    k2 = SS304.thermal_conductivity(T)
    assert np.all(k2 > 0)

def test_material_emissivity():
    assert 0 <= SS310.emissivity <= 1
    assert 0 <= SS304.emissivity <= 1

def test_material_maximum_temperature():
    assert SS310.maximum_temperature > 0
    assert SS304.maximum_temperature > 0

def test_material_youngs_modulus_callable():
    T = np.array([300, 500, 1000])
    E = SS310.youngs_modulus(T)
    assert np.all(E > 0)
    E2 = SS304.youngs_modulus(T)
    assert np.all(E2 > 0)

def test_material_thermal_expansion_coeff_callable():
    T = np.array([300, 500, 1000])
    cte = SS310.thermal_expansion_coeffient(T)
    assert np.all(cte > 0)
    cte2 = SS304.thermal_expansion_coeffient(T)
    assert np.all(cte2 > 0)

def test_material_yield_strength_callable():
    T = np.array([300, 500, 1000])
    ys = SS310.yield_strength(T)
    assert np.all(ys > 0)
    ys2 = SS304.yield_strength(T)
    assert np.all(ys2 > 0)

def test_material_thermal_diffusivity_callable():
    T = np.array([300, 500, 1000])
    alpha = SS310.thermal_diffusivity(T)
    assert np.all(alpha > 0)
    alpha2 = SS304.thermal_diffusivity(T)
    assert np.all(alpha2 > 0)
