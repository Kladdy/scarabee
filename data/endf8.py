import numpy as np
import libraries
from libraries import Nuclide, ThermalScatteringMaterial
import tools.frendy as fdy
import tools.depletion_chain as dc
import h5py


####################################################################################
##                         MODIFY LINES IN THIS BLOCK
library = libraries.endf71_library # Change to desired library
library.download_library_if_not_exists()
# or, define custom library:
# library = libraries.ENDFLibrary(
#     name="Custom Library",
#     label="custom_label",
#     base_path=libraries.base_endf_path / "custom_library",
#     neutrons_path_suffix="neutrons",
#     tsl_path_suffix="thermal_scatt",
# )

temps = [293.0, 500.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0]

hdf5_file_name = f"{library.label}_shem281.h5"

# Set the default group strucutre
fdy.set_default_group_structure("SHEM-281")
fdy.set_default_max_legendre_moments(3)

####################################################################################


# First we initialize the file with the depletion chain data.
# Must be done first as the save method will delete file first if it already exists !
dep_chain = dc.build_depletion_chain("chain_casl_pwr.xml")
dep_chain.save(hdf5_file_name)

# Initialize HDF5 with the library information
print(fdy.get_default_group_structure())
h5 = h5py.File(hdf5_file_name, "r+")
h5.attrs["group-structure"] = fdy.get_default_group_structure().name
h5.attrs["group-bounds"] = fdy.get_default_group_structure().bounds
h5.attrs["ngroups"] = fdy.get_default_group_structure().ngroups
h5.attrs["condensation-scheme"] = fdy.get_default_group_structure().condensation_scheme
h5.attrs["cmfd-condensation-scheme"] = (
    fdy.get_default_group_structure().cmfd_condensation_scheme
)
h5.attrs["first-resonance-group"] = fdy.get_default_group_structure().first_res_grp
h5.attrs["last-resonance-group"] = fdy.get_default_group_structure().last_res_grp
h5.attrs["library"] = library.name

# Define a helper function to process nuclides
def process_nuclide(nuclide: Nuclide, dilutions: list[float] = [1.0e10]):
    N = fdy.FrendyMG()
    N.name = nuclide.label()
    N.endf_file = str(library.neutrons_path / library.get_filename_of_nuclide(nuclide))
    N.label = f"{nuclide.label()} from {library.name}"
    N.temps = temps
    if dilutions:
        N.dilutions = dilutions
    N.process(h5)
    return N

# Define a helper function to process TSL materials
def process_tsl_material(tsl_material: ThermalScatteringMaterial, nuclide: Nuclide, chi: np.ndarray | None, dilutions: list[float] = [1.0e10]):
    tsl_temps = library.get_tsl_temperatures(tsl_material)
    
    N = fdy.FrendyMG()
    N.name = tsl_material.label()
    N.endf_file = str(library.neutrons_path / library.get_filename_of_nuclide(nuclide))
    N.tsl_file = str(library.tsl_path / library.get_filename_of_tsl_material(tsl_material))
    N.tsl_type = tsl_material.subject.lower()
    N.label = f"{tsl_material.label()} from {library.name}"
    N.temps = tsl_temps
    if dilutions:
        N.dilutions = dilutions
    N.process(h5, chi)
    return N

# # We process U235 first so that we can steal its fission spectrum for
# # performing transport correciton calculations
# N = process_nuclide(Nuclide("U", 92, 235))

# if N.chi is not None:
#     chi = N.chi[:]
# else:
#     raise RuntimeError("U235 fission spectrum not found")

chi=None

# Process TSL based evaluations
process_tsl_material(
    ThermalScatteringMaterial("H", "H2O", "hh2o"),
    Nuclide("H", 1, 1),
    chi=chi,
    dilutions=[1.0e10],
)

process_tsl_material(
    ThermalScatteringMaterial("H", "D2O", "dd2o"),
    Nuclide("H", 1, 2),
    chi=chi,
    dilutions=[1.0e10],
)

# Free Gas Nuclides
nuclides = [
    Nuclide("H", 1, 1),
    Nuclide("H", 1, 2),
    Nuclide("He", 2, 3),
    Nuclide("He", 2, 4),
    Nuclide("Li", 3, 6),
    Nuclide("Li", 3, 7),
    Nuclide("Be", 4, 9),
    Nuclide("B", 5, 10),
    Nuclide("B", 5, 11),
    Nuclide("C", 6, 12),
    Nuclide("C", 6, 13),
    Nuclide("N", 7, 14),
    Nuclide("N", 7, 15),
    Nuclide("O", 8, 16),
    Nuclide("O", 8, 17),
    Nuclide("O", 8, 18),
    Nuclide("Na", 11, 23),
    Nuclide("Mg", 12, 24),
    Nuclide("Mg", 12, 25),
    Nuclide("Mg", 12, 26),
    Nuclide("Al", 13, 27),
    Nuclide("Si", 14, 28),
    Nuclide("Si", 14, 29),
    Nuclide("Si", 14, 30),
    Nuclide("Ar", 18, 36),
    Nuclide("Ar", 18, 38),
    Nuclide("Ar", 18, 40),
    Nuclide("Ti", 22, 46),
    Nuclide("Ti", 22, 47),
    Nuclide("Ti", 22, 48),
    Nuclide("Ti", 22, 49),
    Nuclide("Ti", 22, 50),
    Nuclide("Cr", 24, 50),
    Nuclide("Cr", 24, 52),
    Nuclide("Cr", 24, 53),
    Nuclide("Cr", 24, 54),
    Nuclide("Mn", 25, 55),
    Nuclide("Fe", 26, 54),
    Nuclide("Fe", 26, 55),
    Nuclide("Fe", 26, 56),
    Nuclide("Fe", 26, 57),
    Nuclide("Fe", 26, 58),
    Nuclide("Co", 27, 59),
    Nuclide("Ni", 28, 58),
    Nuclide("Ni", 28, 60),
    Nuclide("Ni", 28, 61),
    Nuclide("Ni", 28, 62),
    Nuclide("Ni", 28, 64),
    Nuclide("Cu", 29, 63),
    Nuclide("Cu", 29, 65),
    Nuclide("Zn", 30, 64),
    Nuclide("Zn", 30, 66),
    Nuclide("Zn", 30, 67),
    Nuclide("Zn", 30, 68),
    Nuclide("Zn", 30, 70),
    Nuclide("Br", 35, 81),
    Nuclide("Kr", 36, 82),
    Nuclide("Kr", 36, 83),
    Nuclide("Kr", 36, 84),
    Nuclide("Kr", 36, 85),
    Nuclide("Kr", 36, 86),
    Nuclide("Sr", 38, 89),
    Nuclide("Sr", 38, 90),
    Nuclide("Y", 39, 89),
    Nuclide("Y", 39, 90),
    Nuclide("Y", 39, 91),
    Nuclide("Zr", 40, 90),
    Nuclide("Zr", 40, 91),
    Nuclide("Zr", 40, 92),
    Nuclide("Zr", 40, 93, is_fission_product=True),
    Nuclide("Zr", 40, 94),
    Nuclide("Zr", 40, 95, is_fission_product=True),
    Nuclide("Zr", 40, 96),
    Nuclide("Nb", 41, 95, is_fission_product=True),
    Nuclide("Mo", 42, 92),
    Nuclide("Mo", 42, 94),
    Nuclide("Mo", 42, 95),
    Nuclide("Mo", 42, 96),
    Nuclide("Mo", 42, 97),
    Nuclide("Mo", 42, 98),
    Nuclide("Mo", 42, 99, is_fission_product=True),
    Nuclide("Mo", 42, 100),
    Nuclide("Tc", 43, 99),
    Nuclide("Ru", 44, 99, is_fission_product=True),
    Nuclide("Ru", 44, 100, is_fission_product=True),
    Nuclide("Ru", 44, 101, is_fission_product=True),
    Nuclide("Ru", 44, 102, is_fission_product=True),
    Nuclide("Ru", 44, 103, is_fission_product=True),
    Nuclide("Ru", 44, 104, is_fission_product=True),
    Nuclide("Ru", 44, 105, is_fission_product=True),
    Nuclide("Ru", 44, 106, is_fission_product=True),
    Nuclide("Rh", 45, 103, is_fission_product=True),
    Nuclide("Rh", 45, 105, is_fission_product=True),
    Nuclide("Pd", 46, 104, is_fission_product=True),
    Nuclide("Pd", 46, 105, is_fission_product=True),
    Nuclide("Pd", 46, 106, is_fission_product=True),
    Nuclide("Pd", 46, 107, is_fission_product=True),
    Nuclide("Pd", 46, 108, is_fission_product=True),
    Nuclide("Ag", 47, 107),
    Nuclide("Ag", 47, 109),
    Nuclide("Ag", 47, 110, is_metastable=True),
    Nuclide("Ag", 47, 111),
    Nuclide("Cd", 48, 106),
    Nuclide("Cd", 48, 108),
    Nuclide("Cd", 48, 110),
    Nuclide("Cd", 48, 111),
    Nuclide("Cd", 48, 112),
    Nuclide("Cd", 48, 113),
    Nuclide("Cd", 48, 114),
    Nuclide("Cd", 48, 116),
    Nuclide("In", 49, 113),
    Nuclide("In", 49, 115),
    Nuclide("Sn", 50, 112),
    Nuclide("Sn", 50, 114),
    Nuclide("Sn", 50, 115),
    Nuclide("Sn", 50, 116),
    Nuclide("Sn", 50, 117),
    Nuclide("Sn", 50, 118),
    Nuclide("Sn", 50, 119),
    Nuclide("Sn", 50, 120),
    Nuclide("Sn", 50, 122),
    Nuclide("Sn", 50, 124),
    Nuclide("Sb", 51, 121, is_fission_product=True),
    Nuclide("Sb", 51, 123, is_fission_product=True),
    Nuclide("Sb", 51, 125, is_fission_product=True),
    Nuclide("Te", 52, 127, is_fission_product=True, is_metastable=True),
    Nuclide("Te", 52, 129, is_fission_product=True, is_metastable=True),
    Nuclide("Te", 52, 132, is_fission_product=True),
    Nuclide("I", 53, 127),
    Nuclide("I", 53, 129),
    Nuclide("I", 53, 131, is_fission_product=True),
    Nuclide("I", 53, 135, is_fission_product=True),
    Nuclide("Xe", 54, 128, is_fission_product=True),
    Nuclide("Xe", 54, 129),
    Nuclide("Xe", 54, 130, is_fission_product=True),
    Nuclide("Xe", 54, 131),
    Nuclide("Xe", 54, 132, is_fission_product=True),
    Nuclide("Xe", 54, 133, is_fission_product=True),
    Nuclide("Xe", 54, 134, is_fission_product=True),
    Nuclide("Xe", 54, 135, is_fission_product=True),
    Nuclide("Xe", 54, 136, is_fission_product=True),
    Nuclide("Cs", 55, 133, is_fission_product=True),
    Nuclide("Cs", 55, 134, is_fission_product=True),
    Nuclide("Cs", 55, 135, is_fission_product=True),
    Nuclide("Cs", 55, 136, is_fission_product=True),
    Nuclide("Cs", 55, 137, is_fission_product=True),
    Nuclide("Ba", 56, 134, is_fission_product=True),
    Nuclide("Ba", 56, 137, is_fission_product=True),
    Nuclide("Ba", 56, 140, is_fission_product=True),
    Nuclide("La", 57, 139, is_fission_product=True),
    Nuclide("La", 57, 140, is_fission_product=True),
    Nuclide("Ce", 58, 140, is_fission_product=True),
    Nuclide("Ce", 58, 141, is_fission_product=True),
    Nuclide("Ce", 58, 142, is_fission_product=True),
    Nuclide("Ce", 58, 143, is_fission_product=True),
    Nuclide("Ce", 58, 144, is_fission_product=True),
    Nuclide("Pr", 59, 141, is_fission_product=True),
    Nuclide("Pr", 59, 143, is_fission_product=True),
    Nuclide("Nd", 60, 142, is_fission_product=True),
    Nuclide("Nd", 60, 143, is_fission_product=True),
    Nuclide("Nd", 60, 144, is_fission_product=True),
    Nuclide("Nd", 60, 145, is_fission_product=True),
    Nuclide("Nd", 60, 146, is_fission_product=True),
    Nuclide("Nd", 60, 147, is_fission_product=True),
    Nuclide("Nd", 60, 148, is_fission_product=True),
    Nuclide("Nd", 60, 150, is_fission_product=True),
    Nuclide("Pm", 61, 147, is_fission_product=True),
    Nuclide("Pm", 61, 148, is_fission_product=True),
    Nuclide("Pm", 61, 148, is_fission_product=True, is_metastable=True),
    Nuclide("Pm", 61, 149, is_fission_product=True),
    Nuclide("Pm", 61, 151, is_fission_product=True),
    Nuclide("Sm", 62, 147),
    Nuclide("Sm", 62, 148, is_fission_product=True),
    Nuclide("Sm", 62, 149),
    Nuclide("Sm", 62, 150),
    Nuclide("Sm", 62, 151),
    Nuclide("Sm", 62, 152),
    Nuclide("Sm", 62, 153),
    Nuclide("Sm", 62, 154),
    Nuclide("Eu", 63, 151),
    Nuclide("Eu", 63, 152),
    Nuclide("Eu", 63, 153),
    Nuclide("Eu", 63, 154),
    Nuclide("Eu", 63, 155),
    Nuclide("Eu", 63, 156, is_fission_product=True),
    Nuclide("Gd", 64, 152),
    Nuclide("Gd", 64, 154),
    Nuclide("Gd", 64, 155),
    Nuclide("Gd", 64, 156),
    Nuclide("Gd", 64, 157),
    Nuclide("Gd", 64, 158),
    Nuclide("Gd", 64, 160),
    Nuclide("Tb", 65, 159, is_fission_product=True),
    Nuclide("Tb", 65, 160, is_fission_product=True),
    Nuclide("Tb", 65, 161, is_fission_product=True),
    Nuclide("Dy", 66, 160),
    Nuclide("Dy", 66, 161),
    Nuclide("Dy", 66, 162),
    Nuclide("Dy", 66, 163),
    Nuclide("Dy", 66, 164),
    Nuclide("Ho", 67, 165, is_fission_product=True),
    Nuclide("Er", 68, 162, is_fission_product=True),
    Nuclide("Er", 68, 164, is_fission_product=True),
    Nuclide("Er", 68, 166, is_fission_product=True),
    Nuclide("Er", 68, 167, is_fission_product=True),
    Nuclide("Er", 68, 168, is_fission_product=True),
    Nuclide("Er", 68, 169, is_fission_product=True),
    Nuclide("Er", 68, 170, is_fission_product=True),
    Nuclide("Tm", 69, 169, is_fission_product=True),
    Nuclide("Tm", 69, 170, is_fission_product=True),
    Nuclide("Tm", 69, 171, is_fission_product=True),
    Nuclide("Hf", 72, 174),
    Nuclide("Hf", 72, 176),
    Nuclide("Hf", 72, 177),
    Nuclide("Hf", 72, 178),
    Nuclide("Hf", 72, 179),
    Nuclide("Hf", 72, 180),
    Nuclide("Hf", 72, 181),
    Nuclide("Ta", 73, 181),
    Nuclide("Ta", 73, 182, is_fission_product=True),
    Nuclide("Th", 90, 230),
    Nuclide("Th", 90, 231),
    Nuclide("Th", 90, 232),
    Nuclide("Th", 90, 234),
    Nuclide("Pa", 91, 231),
    Nuclide("Pa", 91, 232),
    Nuclide("Pa", 91, 233),
    Nuclide("U", 92, 232),
    Nuclide("U", 92, 233),
    Nuclide("U", 92, 234),
    Nuclide("U", 92, 236),
    Nuclide("U", 92, 237),
    Nuclide("U", 92, 238),
    Nuclide("Np", 93, 236),
    Nuclide("Np", 93, 237),
    Nuclide("Np", 93, 238),
    Nuclide("Np", 93, 239),
    Nuclide("Pu", 94, 236),
    Nuclide("Pu", 94, 237),
    Nuclide("Pu", 94, 238),
    Nuclide("Pu", 94, 239),
    Nuclide("Pu", 94, 240),
    Nuclide("Pu", 94, 241),
    Nuclide("Pu", 94, 242),
    Nuclide("Pu", 94, 243),
    Nuclide("Pu", 94, 244),
    Nuclide("Am", 95, 241),
    Nuclide("Am", 95, 242, is_metastable=True),
    Nuclide("Am", 95, 243),
    Nuclide("Cm", 96, 242),
    Nuclide("Cm", 96, 243),
    Nuclide("Cm", 96, 244),
    Nuclide("Cm", 96, 245),
    Nuclide("Cm", 96, 246),
]

[process_nuclide(nuclide) for nuclide in nuclides]

h5.close()
