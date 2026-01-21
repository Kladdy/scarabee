from scarabee import (
    NDLibrary,
    MaterialComposition,
    Material,
    Fraction,
    DensityUnits,
    set_output_file,
)
from scarabee.reseau import FuelPin, GuideTube, PWRAssembly, Reflector
from scarabee.coeur import SimpleTile, QuadrantsTile
import pickle

name = "F31_0"

set_output_file(name + "_out.txt")

ndl = NDLibrary()

L = 3

# Define all Materials
Fuel31Comp = MaterialComposition(Fraction.Atoms, name="Fuel 3.1%")
Fuel31Comp.add_leu(3.1, 1.0)
Fuel31Comp.add_element("O", 2.0)
Fuel31 = Material(Fuel31Comp, 575.0, 10.30166, DensityUnits.g_cm3, ndl)
Fuel31.max_legendre_order = L

CladComp = MaterialComposition(Fraction.Weight, name="Zircaloy 4")
CladComp.add_element("O", 0.00125)
CladComp.add_element("Cr", 0.0010)
CladComp.add_element("Fe", 0.0021)
CladComp.add_element("Zr", 0.98115)
CladComp.add_element("Sn", 0.0145)
Clad = Material(CladComp, 575.0, 6.55, DensityUnits.g_cm3, ndl)
Clad.max_legendre_order = L

HeComp = MaterialComposition(Fraction.Atoms, name="He Gas")
HeComp.add_element("He", 1.0)
He = Material(HeComp, 575.0, 0.0015981, DensityUnits.g_cm3, ndl)
He.max_legendre_order = L

# Define a guide tube
gt = GuideTube(inner_radius=0.56134, outer_radius=0.60198, clad=Clad)

# Define fuel pin
fp = FuelPin(
    fuel=Fuel31,
    fuel_radius=0.39218,
    gap=He,
    gap_radius=0.40005,
    clad=Clad,
    clad_radius=0.45720,
)

cells = [
    [fp, fp, fp, fp, fp, fp, fp, fp, fp],
    [fp, fp, fp, fp, fp, fp, fp, fp, fp],
    [gt, fp, fp, gt, fp, fp, fp, fp, fp],
    [fp, fp, fp, fp, fp, gt, fp, fp, fp],
    [fp, fp, fp, fp, fp, fp, fp, fp, fp],
    [gt, fp, fp, gt, fp, fp, gt, fp, fp],
    [fp, fp, fp, fp, fp, fp, fp, fp, fp],
    [fp, fp, fp, fp, fp, fp, fp, fp, fp],
    [gt, fp, fp, gt, fp, fp, gt, fp, fp],
]

# Define assembly
asmbly = PWRAssembly(
    pitch=1.25984,
    assembly_pitch=21.50364,
    shape=(17, 17),
    moderator={'boron-ppm': 975., 'temperature': 575., 'pressure': 15.5132, 'legendre-order': L},
    independent_quadrant=True,
    cells=cells,
    ndl=ndl,
)
asmbly.solve()

ct = QuadrantsTile.from_independent_quadrant(asmbly.diffusion_data, asmbly.form_factors)
pickle.dump(ct, open(f'{name}.pkl', 'wb'))

# Homogenize the assembly
avg_fuel_asmbly = asmbly._asmbly_moc.homogenize()

# Baffle material
SS304Comp = MaterialComposition(Fraction.Weight, name="SS304")
SS304Comp.add_element("Si", 0.0060)
SS304Comp.add_element("Cr", 0.1900)
SS304Comp.add_element("Mn", 0.0200)
SS304Comp.add_element("Fe", 0.6840)
SS304Comp.add_element("Ni", 0.1000)
SS304 = Material(SS304Comp, 575.0, 8.03, DensityUnits.g_cm3, ndl)
SS304.max_legendre_order = 3

refl = Reflector(
    avg_fuel_asmbly,
    moderator=asmbly._moderator_xs,
    assembly_width=21.50364,
    gap_width=0.1627,
    baffle_width=2.2225,
    baffle=SS304,
    ndl=ndl,
)
refl.anisotropic = True
refl.solve()

rt = SimpleTile(refl.diffusion_data, refl.form_factors)
pickle.dump(rt, open('reflector.pkl', 'wb'))
